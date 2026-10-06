// fitDecay.cpp
// Usage: ./fitDecay <file.root>
// Averages the neutron-like waveforms that hold exactly one pulse, and fits the
// average to extract the physical shape of the pulse: its decay times and
// weights, and its rise time. It is done for all selected events together and
// for each bin of pulse height (BIN_EDGES), to show if the decay depends on the
// size of the pulse.
//
// Selection: max amplitude between the first and last of BIN_EDGES, neutron-like
// (at least NEUTRON_MIN_STABLE_SAMPLES samples in a row above
// NEUTRON_AMPLITUDE_THRESHOLD), PSD ratio in [PSD_MIN, PSD_MAX], exactly one
// pulse (see PileupFinder.h), and that pulse starting where the digitizer
// triggered (between ONSET_MIN_SAMPLE and ONSET_MAX_SAMPLE). The last condition
// leaves out events in which a gamma spike later in the record was found as the
// only pulse while the neutron that triggered the digitizer was too small to be
// found.
//
// Each average is fitted with three exponential components, from its peak
// (plus FIT_START_AFTER_PEAK_NS) to FIT_END_NS, and the decay times and weights
// are printed. One PDF per average is written: the average pulse (peak
// normalised to 1) with the fit as a dashed line.

#include "Waveform.h"
#include "PSD.h"
#include "PileupFinder.h"
#include "DecayAnalysis.h"
#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>

static const double NS_PER_SAMPLE = 2.0;

// Pulse-height bins (max amplitude, ADC). The first and last edge are the
// selection window.
static const double BIN_EDGES[] = {4000.0, 6000.0, 8000.0, 12000.0};

// Neutron-like: at least NEUTRON_MIN_STABLE_SAMPLES samples in a row above
// NEUTRON_AMPLITUDE_THRESHOLD.
static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;

// PSD gates, same as runPSDSamples13, and the PSD cut.
static const double POST_TRIGGER_PERCENT = 80.0;
static const double N_COEFFICIENT        = 8.0;
static const double CONSTANT_LATENCY     = 73.2889;
static const double SHORT_GATE_NS        = 40.0;
static const double LONG_GATE_NS         = 300.0;
static const double START_SHIFT_NS       = 7.0;
static const double PSD_MIN = 0.7;
static const double PSD_MAX = 1.0;

// Pulse finder, same as findPileup.
static const double CLIMB_ADC      = 500.0;
static const double MIN_FRACTION   = 0.3;
static const double NOISE_SIGMAS   = 8.0;
static const int    NOISE_SAMPLES  = 30;
static const int    WINDOW_SAMPLES = 10;
static const int    PEAK_SAMPLES   = 25;
static const int    BLANK_SAMPLES  = 50;
static const int    MAX_PULSES     = 8;

// The pulse must start where the digitizer triggered (samples, as found by the
// pulse finder; 310-360 is 620-720 ns).
static const int ONSET_MIN_SAMPLE = 310;
static const int ONSET_MAX_SAMPLE = 360;

// Averaging and fitting, see DecayAnalysis.h.
static const double ALIGN_FRACTION    = 0.5;
static const int    ALIGN_SEARCH_BEFORE = 2;
static const int    ALIGN_SEARCH_AFTER  = 40;
static const double AREA_NS           = 1500.0;
static const double PRE_NS            = 200.0;
static const double POST_NS           = 3300.0;
static const double FIT_END_NS        = 3300.0;
static const double FIT_START_AFTER_PEAK_NS = 0.0;
static const double SYSTEMATIC_FLOOR  = 0.001;
static const double SIGMA_STARTS[]    = {2.0, 6.0};
static const double FAST_STARTS[]     = {30.0, 80.0};
static const double MEDIUM_STARTS[]   = {300.0, 800.0};
static const double SLOW_STARTS[]     = {2000.0, 8000.0};
static const int    MAX_EVALUATIONS   = 3000;
static const double TOLERANCE         = 1e-10;
static const long   MIN_EVENTS        = 200;

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    std::string inFile = Utilities::INPUT_DIR + argv[1];

    const int nBins = (int)(sizeof(BIN_EDGES) / sizeof(BIN_EDGES[0])) - 1;
    const double low  = BIN_EDGES[0];
    const double high = BIN_EDGES[nBins];

    PileupFinder::Settings fs;
    fs.jump          = CLIMB_ADC;
    fs.minFraction   = MIN_FRACTION;
    fs.noiseSigmas   = NOISE_SIGMAS;
    fs.noiseSamples  = NOISE_SAMPLES;
    fs.windowSamples = WINDOW_SAMPLES;
    fs.peakSamples   = PEAK_SAMPLES;
    fs.blankSamples  = BLANK_SAMPLES;
    fs.maxPulses     = MAX_PULSES;

    DecayAnalysis::Settings ds;
    ds.nsPerSample     = NS_PER_SAMPLE;
    ds.alignFraction   = ALIGN_FRACTION;
    ds.searchBefore    = ALIGN_SEARCH_BEFORE;
    ds.searchAfter     = ALIGN_SEARCH_AFTER;
    ds.peakSamples     = PEAK_SAMPLES;
    ds.areaNs          = AREA_NS;
    ds.preNs           = PRE_NS;
    ds.postNs          = POST_NS;
    ds.fitEndNs        = FIT_END_NS;
    ds.systematicFloor = SYSTEMATIC_FLOOR;
    for (double v : SIGMA_STARTS)  ds.sigmaStarts.push_back(v);
    for (double v : FAST_STARTS)   ds.fastStarts.push_back(v);
    for (double v : MEDIUM_STARTS) ds.mediumStarts.push_back(v);
    for (double v : SLOW_STARTS)   ds.slowStarts.push_back(v);
    ds.maxEvaluations  = MAX_EVALUATIONS;
    ds.tolerance       = TOLERANCE;
    ds.minEvents       = MIN_EVENTS;

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
    tree->SetBranchAddress("recordLength", &recordLength);
    tree->SetBranchAddress("waveform", waveform);

    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "No entries in file: " << inFile << std::endl;
        return 1;
    }

    PileupFinder finder(fs);
    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, 0.0, START_SHIFT_NS);

    // Set 0 holds all events; set 1 + b holds the events in bin b.
    std::vector<DecayAnalysis> sets;
    std::vector<std::string> labels;
    sets.push_back(DecayAnalysis(ds));
    labels.push_back("all");
    for (int b = 0; b < nBins; b++) {
        char label[64];
        snprintf(label, sizeof(label), "%.0f-%.0f", BIN_EDGES[b], BIN_EDGES[b + 1]);
        sets.push_back(DecayAnalysis(ds));
        labels.push_back(label);
    }

    std::string outDir = Utilities::makeOutputDir(argv[0], argv[1]);
    gSystem->mkdir(outDir.c_str(), true);

    long nWindow = 0, nPsdCut = 0, nNotSingle = 0, nOffTrigger = 0, nNoAlign = 0, nUsed = 0;

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
        nWindow++;

        const std::vector<double>& shape = wFull.getSubtractedADC();
        PSD::Result r = psd.analyzeAt(shape, recordLength, startSample);
        if (!(r.psdRatio >= PSD_MIN && r.psdRatio <= PSD_MAX)) { nPsdCut++; continue; }

        std::vector<PileupFinder::Pulse> pulses = finder.find(shape);
        if (pulses.size() != 1) { nNotSingle++; continue; }

        int onset = pulses[0].onsetSample;
        if (onset < ONSET_MIN_SAMPLE || onset > ONSET_MAX_SAMPLE) { nOffTrigger++; continue; }

        if (!sets[0].add(shape, onset)) { nNoAlign++; continue; }
        nUsed++;

        for (int b = 0; b < nBins; b++) {
            bool inBin = (amp >= BIN_EDGES[b]) && (amp < BIN_EDGES[b + 1] || (b == nBins - 1 && amp <= BIN_EDGES[b + 1]));
            if (inBin) sets[1 + b].add(shape, onset);
        }
    }

    std::cout << inFile << ": " << nEntries << " events" << std::endl;
    std::cout << "Neutron-like, max amplitude in [" << low << ", " << high << "] ADC: " << nWindow << std::endl;
    std::cout << "  removed by the PSD cut [" << PSD_MIN << ", " << PSD_MAX << "]:       " << nPsdCut << std::endl;
    std::cout << "  removed, not exactly one pulse:               " << nNotSingle << std::endl;
    std::cout << "  removed, pulse not at the trigger position:   " << nOffTrigger << std::endl;
    std::cout << "  removed, could not be lined up:               " << nNoAlign << std::endl;
    std::cout << "  used:                                         " << nUsed << std::endl;

    for (size_t s = 0; s < sets.size(); s++) {
        if (sets[s].numEvents() < MIN_EVENTS) {
            std::cout << "\n" << labels[s] << " ADC: " << sets[s].numEvents()
                      << " waveforms, too few to fit (fewer than " << MIN_EVENTS << ")" << std::endl;
            continue;
        }

        DecayAnalysis::FitResult fit = sets[s].fit(3, false, sets[s].peakTimeNs() + FIT_START_AFTER_PEAK_NS, true);
        DecayAnalysis::print(fit, labels[s] + " ADC");

        std::string pdf = outDir + "decay_" + labels[s] + ".pdf";
        sets[s].writePdf(pdf, fit);
        std::cout << "  plot: " << pdf << std::endl;
    }

    f->Close();
    return 0;
}
