// findPileup.cpp
// Usage: ./findPileup <file.root> [low_ADC] [high_ADC]
// Takes the neutron-like events whose max amplitude is in [low, high] (default
// 4000-12000 ADC) and whose PSD ratio is in [PSD_MIN, PSD_MAX], finds every
// pulse in each one, and prints how many events hold 1, 2, 3 or more pulses.
//
// A climb is a place where the average of the next 20 ns is at least CLIMB_ADC
// above the average of the previous 20 ns. The height of a climb is its peak
// (within 50 ns) minus the level just before it. The first climb is always a
// pulse (it triggered the digitizer). A later climb is a pulse if its height is
// at least MIN_FRACTION of the event's biggest peak AND at least NOISE_SIGMAS
// times the noise just before it (the scatter of the signal around a straight
// line over the previous 60 ns), so the limit follows the size and the noise of
// each event. Climbs within 100 ns of a pulse are ignored, because the top of a
// pulse is jagged.
//
// Every event with 2 or more pulses is written to PDFs with the pulses marked.
// So is every GOOD_EVERY-th event with exactly 1 pulse (up to GOOD_MAX of them),
// spread over the whole run, to show what a clean event looks like.
// A second, wide-open search records every climb after the first pulse
// (pulseHeights.pdf), and five statistics plots (pileupStatistics.pdf) show
// whether second pulses are random coincidences or tied to the first pulse.

#include "Waveform.h"
#include "PSD.h"
#include "PileupFinder.h"
#include "PileupUtilities.h"
#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>

static const double DEFAULT_LOW_ADC  = 4000.0;
static const double DEFAULT_HIGH_ADC = 12000.0;

static const double NS_PER_SAMPLE = 2.0;

// Neutron-like: at least NEUTRON_MIN_STABLE_SAMPLES samples in a row above
// NEUTRON_AMPLITUDE_THRESHOLD.
static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;

// PSD gates, same as runPSDSamples13. The PSD ratio is used for the cut below
// and for the label on each plot.
static const double POST_TRIGGER_PERCENT = 80.0;
static const double N_COEFFICIENT        = 8.0;
static const double CONSTANT_LATENCY     = 73.2889;
static const double SHORT_GATE_NS        = 40.0;
static const double LONG_GATE_NS         = 300.0;
static const double START_SHIFT_NS       = 7.0;

// Events whose PSD ratio is outside [PSD_MIN, PSD_MAX] are left out.
static const double PSD_MIN = 0.7;
static const double PSD_MAX = 1.0;

// Pulse finder, see PileupFinder.h.
static const double CLIMB_ADC      = 500.0;    // smallest climb to look at (the hardware trigger level)
static const double MIN_FRACTION   = 0.3;      // a later pulse is at least this fraction of the event's peak ...
static const double NOISE_SIGMAS   = 8.0;      // ... and at least this many times the noise before it
static const int    NOISE_SAMPLES  = 30;       // the noise is measured over the 60 ns before the climb
static const int    WINDOW_SAMPLES = 10;       // 20 ns on each side
static const int    PEAK_SAMPLES   = 25;       // 50 ns
static const int    BLANK_SAMPLES  = 50;       // 100 ns
static const int    MAX_PULSES     = 8;

// Wide-open finder for the height study: every climb counts, so every climb is
// seen. Climbs within GUARD_SAMPLES of the first pulse are not counted in the
// study, because the top of a pulse is jagged.
static const int    WIDE_BLANK_SAMPLES = 10;   // 20 ns
static const int    WIDE_MAX_PULSES    = 60;
static const int    GUARD_SAMPLES      = 50;   // 100 ns

// Statistics plots (pileupStatistics.pdf), see PileupUtilities.h.
static const int    FIRST_BINS   = 60;
static const double FIRST_MAX    = 12000.0;
static const int    LATER_BINS   = 60;
static const double LATER_MAX    = 6000.0;
static const int    GAP_BINS     = 68;
static const double GAP_MAX      = 4080.0;
static const double GAP_FLAT_END = 3400.0;
static const int    FRAC_BINS    = 60;
static const double FRAC_MAX     = 1.2;
static const int    SIGMA_BINS   = 60;
static const double SIGMA_MAX    = 60.0;

// Events with exactly 1 pulse that are written to PDFs: every GOOD_EVERY-th one,
// up to GOOD_MAX of them.
static const long   GOOD_EVERY = 50;
static const long   GOOD_MAX   = 3000;

// Plots and PDFs, see PileupUtilities.h.
static const int    PAGES_PER_PDF     = 500;
static const long   MAX_STORED        = 20000;
static const int    CANVAS_WIDTH      = 1000;
static const int    CANVAS_HEIGHT     = 600;
static const double HEADROOM          = 1.35;
static const int    MAX_PULSES_LISTED = 8;
static const int    HEIGHT_BINS       = 120;
static const double HEIGHT_MAX        = 6000.0;
static const int    RATIO_BINS        = 75;
static const double RATIO_MAX         = 1.5;
static const double TABLE_EDGES[]     = {1000.0, 1500.0, 2000.0, 3000.0, 4000.0};

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;

    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> [low_ADC] [high_ADC]" << std::endl;
        return 1;
    }

    std::string inFile = Utilities::INPUT_DIR + argv[1];
    double low  = (argc >= 3) ? std::atof(argv[2]) : DEFAULT_LOW_ADC;
    double high = (argc >= 4) ? std::atof(argv[3]) : DEFAULT_HIGH_ADC;

    if (low >= high) {
        std::cerr << "low (" << low << ") must be smaller than high (" << high << ")" << std::endl;
        return 1;
    }

    PileupFinder::Settings fs;
    fs.jump          = CLIMB_ADC;
    fs.minFraction   = MIN_FRACTION;
    fs.noiseSigmas   = NOISE_SIGMAS;
    fs.noiseSamples  = NOISE_SAMPLES;
    fs.windowSamples = WINDOW_SAMPLES;
    fs.peakSamples   = PEAK_SAMPLES;
    fs.blankSamples  = BLANK_SAMPLES;
    fs.maxPulses     = MAX_PULSES;

    PileupFinder::Settings ws = fs;
    ws.minFraction   = 0.0;
    ws.noiseSigmas   = 0.0;
    ws.blankSamples  = WIDE_BLANK_SAMPLES;
    ws.maxPulses     = WIDE_MAX_PULSES;

    PileupUtilities::Settings ps;
    ps.nsPerSample     = NS_PER_SAMPLE;
    ps.pagesPerPdf     = PAGES_PER_PDF;
    ps.maxStored       = MAX_STORED;
    ps.canvasWidth     = CANVAS_WIDTH;
    ps.canvasHeight    = CANVAS_HEIGHT;
    ps.headroom        = HEADROOM;
    ps.maxPulsesListed = MAX_PULSES_LISTED;
    ps.heightBins      = HEIGHT_BINS;
    ps.heightMax       = HEIGHT_MAX;
    ps.ratioBins       = RATIO_BINS;
    ps.ratioMax        = RATIO_MAX;
    ps.fracBins        = FRAC_BINS;
    ps.fracMax         = FRAC_MAX;
    ps.sigmaBins       = SIGMA_BINS;
    ps.sigmaMax        = SIGMA_MAX;
    ps.cutFraction     = MIN_FRACTION;
    ps.cutSigmas       = NOISE_SIGMAS;
    ps.firstBins       = FIRST_BINS;
    ps.firstMax        = FIRST_MAX;
    ps.laterBins       = LATER_BINS;
    ps.laterMax        = LATER_MAX;
    ps.gapBins         = GAP_BINS;
    ps.gapMax          = GAP_MAX;
    ps.gapFlatEnd      = GAP_FLAT_END;
    for (double e : TABLE_EDGES) ps.tableEdges.push_back(e);

    TFile* f = TFile::Open(inFile.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open file: " << inFile << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)f->Get("waveforms");
    if (!tree) {
        std::cerr << "Could not find tree 'waveforms' in file: " << inFile << std::endl;
        return 1;
    }

    Int_t recordLength;
    int waveform[20000];
    Double_t unixTime;
    tree->SetBranchAddress("recordLength", &recordLength);
    tree->SetBranchAddress("waveform", waveform);
    tree->SetBranchAddress("unixTime", &unixTime);

    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "No entries in file: " << inFile << std::endl;
        return 1;
    }

    tree->GetEntry(0);
    double tStart = unixTime;
    tree->GetEntry(nEntries - 1);
    double duration = unixTime - tStart;

    PileupFinder finder(fs);
    PileupFinder wide(ws);
    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, 0.0, START_SHIFT_NS);

    std::string outDir = Utilities::makeOutputDir(argv[0], argv[1]);
    gSystem->mkdir(outDir.c_str(), true);
    PileupUtilities plots(outDir, ps);

    long nSelected = 0, nNone = 0, n1 = 0, n2 = 0, n3 = 0, n4plus = 0, totalPulses = 0;
    long nGoodPlotted = 0, nPsdCut = 0;
    std::vector<double> gaps;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);

        Waveform w(waveform, recordLength, false);
        double amp = w.getMaxAmpADC();
        if (amp < low || amp > high) continue;

        Waveform wFull(waveform, recordLength, true);
        int startSample;
        double localPeak;
        if (!wFull.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                           startSample, localPeak)) continue;

        const std::vector<double>& shape = wFull.getSubtractedADC();
        PSD::Result r = psd.analyzeAt(shape, recordLength, startSample);
        if (!(r.psdRatio >= PSD_MIN && r.psdRatio <= PSD_MAX)) {
            nPsdCut++;
            continue;
        }

        nSelected++;
        std::vector<PileupFinder::Pulse> pulses = finder.find(shape);
        plots.addEvent(wide.find(shape), pulses, amp, GUARD_SAMPLES);

        int np = (int)pulses.size();
        totalPulses += np;
        if      (np == 0) nNone++;
        else if (np == 1) n1++;
        else if (np == 2) n2++;
        else if (np == 3) n3++;
        else              n4plus++;

        for (int k = 1; k < np; k++) {
            gaps.push_back((pulses[k].onsetSample - pulses[k - 1].onsetSample) * NS_PER_SAMPLE);
        }

        bool plotThis = (np >= 2);
        if (np == 1 && n1 % GOOD_EVERY == 0 && nGoodPlotted < GOOD_MAX) {
            plotThis = true;
            nGoodPlotted++;
        }
        if (plotThis) plots.add(i, shape, pulses, r.psdRatio);
    }

    std::cout << inFile << ": " << nEntries << " total events, run duration = " << duration << " s" << std::endl;
    std::cout << "Selection: neutron-like, PSD in [" << PSD_MIN << ", " << PSD_MAX
              << "], max amplitude in [" << low << ", " << high << "] ADC" << std::endl;
    std::cout << "Finder: a climb is where the average of the next " << WINDOW_SAMPLES * NS_PER_SAMPLE
              << " ns is at least " << CLIMB_ADC << " ADC above the average of the previous "
              << WINDOW_SAMPLES * NS_PER_SAMPLE << " ns;" << std::endl;
    std::cout << "        the first climb is always a pulse; a later climb needs a height of at least "
              << MIN_FRACTION << " of the event's peak and at least " << NOISE_SIGMAS
              << " times the noise before it (measured over " << NOISE_SAMPLES * NS_PER_SAMPLE << " ns);" << std::endl;
    std::cout << "        climbs within " << BLANK_SAMPLES * NS_PER_SAMPLE << " ns of a pulse are ignored" << std::endl;
    std::cout << "\nRemoved by the PSD cut: " << nPsdCut << std::endl;
    std::cout << "Selected events: " << nSelected;
    if (duration > 0) std::cout << "  (" << nSelected / duration << " Hz)";
    std::cout << std::endl;

    std::cout << "Pulses found per event:" << std::endl;
    std::cout << "  0 pulses:  " << nNone << std::endl;
    std::cout << "  1 pulse:   " << n1 << std::endl;
    std::cout << "  2 pulses:  " << n2 << "  (" << (nSelected > 0 ? 100.0 * n2 / nSelected : 0.0) << "%)" << std::endl;
    std::cout << "  3 pulses:  " << n3 << "  (" << (nSelected > 0 ? 100.0 * n3 / nSelected : 0.0) << "%)" << std::endl;
    std::cout << "  4+ pulses: " << n4plus << "  (" << (nSelected > 0 ? 100.0 * n4plus / nSelected : 0.0) << "%)" << std::endl;

    long nMulti = n2 + n3 + n4plus;
    std::cout << "Events with 2 or more pulses: " << nMulti << "  ("
              << (nSelected > 0 ? 100.0 * nMulti / nSelected : 0.0) << "%)" << std::endl;
    std::cout << "Total pulses found: " << totalPulses;
    if (duration > 0) std::cout << "  (" << totalPulses / duration << " Hz, against "
                                << nSelected / duration << " Hz of events)";
    std::cout << std::endl;

    if (!gaps.empty()) {
        std::sort(gaps.begin(), gaps.end());
        std::cout << "Median gap between consecutive pulses: " << gaps[gaps.size() / 2] << " ns" << std::endl;
    }

    std::cout << std::endl;
    plots.printHeightSummary();
    plots.writeHeightPlots(outDir + "pulseHeights.pdf");
    std::cout << "Height plots: " << outDir << "pulseHeights.pdf" << std::endl;
    plots.writeStatisticsPlots(outDir + "pileupStatistics.pdf");
    std::cout << "Statistics plots: " << outDir << "pileupStatistics.pdf" << std::endl;

    std::cout << std::endl;
    plots.write();
    if (plots.numDropped() > 0) {
        std::cout << "Note: " << plots.numDropped() << " events were counted but not plotted "
                  << "(MAX_STORED = " << MAX_STORED << ")" << std::endl;
    }

    f->Close();
    return 0;
}
