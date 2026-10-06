// countAmplitudeWindow.cpp
// Usage: ./countAmplitudeWindow <file.root> [low_ADC] [high_ADC]
// Counts events whose max amplitude falls in [low, high] (default 4000-12000
// ADC) and reports the rate with its statistical error. Events in the window
// are also split, using the duration tests from runPSDSamples13, into
// neutron-like (long stretch above threshold), single pulse, and pile-up.
// The neutron-like events are further split by whether their PSD ratio, with
// gates anchored at the detected pulse start, is inside [PSD_MIN, PSD_MAX].
// Console output only.

#include "Waveform.h"
#include "PSD.h"
#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>

static const double DEFAULT_LOW_ADC  = 4000.0;
static const double DEFAULT_HIGH_ADC = 12000.0;

static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;
static const int    PILEUP_MIN_QUIET_SAMPLES    = 15;

static const double PSD_MIN = 0.7;
static const double PSD_MAX = 1.0;

// PSD gate settings, same as runPSDSamples13.
static const double POST_TRIGGER_PERCENT = 80.0;
static const double N_COEFFICIENT        = 8.0;
static const double CONSTANT_LATENCY     = 73.2889;
static const double NS_PER_SAMPLE        = 2.0;
static const double SHORT_GATE_NS        = 40.0;
static const double LONG_GATE_NS         = 300.0;
static const double START_SHIFT_NS       = 7.0;

static void printClass(const char* label, long n, long nInWindow, double duration)
{
    double frac = (nInWindow > 0) ? 100.0 * n / nInWindow : 0.0;
    std::cout << "  " << label << n << " (" << frac << "% of in-window)";
    if (duration > 0) {
        std::cout << "  " << n / duration << " +/- " << std::sqrt((double)n) / duration << " Hz";
    }
    std::cout << std::endl;
}

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
    double tEnd = unixTime;
    double duration = tEnd - tStart;

    long nBelow = 0, nInWindow = 0, nAbove = 0;
    long nNeutronLike = 0, nSingle = 0, nPileup = 0;
    long nNeutronPsdIn = 0, nNeutronPsdOut = 0;

    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, 0.0, START_SHIFT_NS);

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, false);
        double amp = w.getMaxAmpADC();

        if (amp < low) {
            nBelow++;
        } else if (amp <= high) {
            nInWindow++;

            // Full waveform is needed for the classification, so build it
            // only for the (few) events inside the window.
            Waveform wFull(waveform, recordLength, true);
            int startSample;
            double localPeak;
            if (wFull.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                              startSample, localPeak)) {
                nNeutronLike++;
                PSD::Result r = psd.analyzeAt(wFull.getSubtractedADC(), recordLength, startSample);
                if (r.psdRatio >= PSD_MIN && r.psdRatio <= PSD_MAX) nNeutronPsdIn++;
                else                                                nNeutronPsdOut++;
            } else {
                auto pulses = wFull.countSeparatePulses(NEUTRON_AMPLITUDE_THRESHOLD, PILEUP_MIN_QUIET_SAMPLES);
                if (pulses.size() == 1) nSingle++;
                else                    nPileup++;
            }
        } else {
            nAbove++;
        }
    }

    std::cout << inFile << ": " << nEntries << " total events, run duration = "
              << duration << " s" << std::endl;
    std::cout << "Window: [" << low << ", " << high << "] ADC" << std::endl;
    std::cout << "Below window: " << nBelow << std::endl;
    std::cout << "In window:    " << nInWindow << std::endl;
    std::cout << "Above window: " << nAbove << std::endl;

    if (duration > 0) {
        double rate = nInWindow / duration;
        double err  = std::sqrt((double)nInWindow) / duration;
        std::cout << "Window rate = " << rate << " +/- " << err << " Hz (statistical)" << std::endl;
    } else {
        std::cout << "Run duration is not positive (unixTime missing?), rate not computed" << std::endl;
    }

    std::cout << "\nIn-window events by type:" << std::endl;
    printClass("Neutron-like:  ", nNeutronLike, nInWindow, duration);
    std::cout << "    of which, PSD in [" << PSD_MIN << ", " << PSD_MAX << "]:" << std::endl;
    printClass("  PSD inside:  ", nNeutronPsdIn,  nInWindow, duration);
    printClass("  PSD outside: ", nNeutronPsdOut, nInWindow, duration);
    printClass("Single pulse:  ", nSingle,      nInWindow, duration);
    printClass("Pile-up:       ", nPileup,      nInWindow, duration);

    f->Close();
    return 0;
}