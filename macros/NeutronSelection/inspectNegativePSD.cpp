// inspectNegativePSD.cpp
// Usage: ./inspectNegativePSD <file.root>
//
// Classifies events the same way as runPSDSamples13.cpp (neutron-like /
// gamma-like / pile-up-excluded), computes PSD using the corrected
// detected-pulse anchor (analyzeAt), and prints + plots every event whose
// PSD ratio falls below PSD_MAX_NEGATIVE (0 by default -- i.e. all
// negative-PSD events), up to MAX_EVENTS.

#include "Waveform.h"
#include "PSD.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>

static const std::string OUTPUT_PATH = "results/";
static const std::string INPUT_DIR   = "/Users/david/DTGAnalysis/data/testruns/";

static const double POST_TRIGGER_PERCENT = 80.0;
static const double N_COEFFICIENT        = 8.0;
static const double CONSTANT_LATENCY     = 73.2889;
static const double NS_PER_SAMPLE        = 2.0;

static const double SHORT_GATE_NS  = 40.0;
static const double LONG_GATE_NS   = 300.0;
static const double START_SHIFT_NS = 7.0;

static const double AMPLITUDE_CUTOFF = 50.0;

static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;
static const int    PILEUP_MIN_QUIET_SAMPLES    = 15;

// Only events with PSD ratio below this get printed/plotted (0 = all negatives).
static const double PSD_MAX_NEGATIVE = 0.0;

static const int MAX_EVENTS = 300;

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    std::string inFile = INPUT_DIR + argv[1];

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

    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, 0.0, START_SHIFT_NS);

    std::string outName = argv[1];
    for (char& ch : outName) if (ch == '.') ch = '_';
    std::string pdfName = OUTPUT_PATH + "negativePSD_" + outName + ".pdf";

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetGrid();
    c->Print((pdfName + "[").c_str());

    std::cout << "event\tclass\tmaxAmp\tlocalPeak\tstartSample\tqShort\t\tqLong\t\tPSD" << std::endl;

    Long64_t nEntries = tree->GetEntries();
    int found = 0;

    for (Long64_t i = 0; i < nEntries && found < MAX_EVENTS; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);

        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;

        const std::vector<double>& shape = w.getSubtractedADC();

        int startSample;
        double localPeak;
        std::string cls;

        bool isNeutron = w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                                    startSample, localPeak);
        if (isNeutron) {
            cls = "neutron";
        } else {
            std::vector<int> pulseStarts;
            auto pulses = w.countSeparatePulses(NEUTRON_AMPLITUDE_THRESHOLD, PILEUP_MIN_QUIET_SAMPLES, pulseStarts);
            if (pulses.size() != 1) continue; // pile-up, not part of either sample
            startSample = pulseStarts[0];
            localPeak = pulses[0];
            cls = "gamma";
        }

        PSD::Result r = psd.analyzeAt(shape, recordLength, startSample);
        if (r.psdRatio >= PSD_MAX_NEGATIVE) continue;

        std::cout << i << "\t" << cls << "\t" << w.getMaxAmpADC()
                  << "\t" << localPeak << "\t" << startSample
                  << "\t" << r.qShort << "\t" << r.qLong << "\t" << r.psdRatio << std::endl;

        psd.drawEventWaveformAt(shape, recordLength, startSample, i, c, pdfName);
        found++;
    }

    c->Print((pdfName + "]").c_str());

    std::cout << "\nFound " << found << " negative-PSD events (capped at " << MAX_EVENTS << ")" << std::endl;
    std::cout << "Saved: " << pdfName << std::endl;

    f->Close();
    return 0;
}