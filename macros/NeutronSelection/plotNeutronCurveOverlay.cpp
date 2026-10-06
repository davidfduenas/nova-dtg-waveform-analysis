// plotNeutronCurveOverlay.cpp
// Usage: ./plotNeutronCurveOverlay <file.root> [event1] [event2] ...
// Overlays the paper's reference neutron decay curve, anchored at each
// event's own detected pulse start.
//
// With no event indices given: scans the whole file and overlays the
// curve on the first MAX_SAMPLE_PLOTS neutron-like events found.
// With specific indices given: overlays the curve on exactly those
// events (each must itself pass isNeutronLike(), or it is skipped).

#include "Waveform.h"
#include "PSD.h"
#include "NeutronCurveOverlay.h"
#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>

static const double SAMPLE_SPACING_NS = 2.0;
static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;

static const double POST_TRIGGER_PERCENT = 80.0;
static const double N_COEFFICIENT        = 8.0;
static const double CONSTANT_LATENCY     = 73.2889;
static const double SHORT_GATE_NS        = 40.0;
static const double LONG_GATE_NS         = 300.0;
static const double START_SHIFT_NS       = 7.0;

static const int MAX_SAMPLE_PLOTS = 200;

static void plotOneEvent(TTree* tree, Long64_t idx, int recordLength, int* waveform,
                          PSD& psd, NeutronCurveOverlay& overlay,
                          TCanvas* c, const std::string& pdfName)
{
    tree->GetEntry(idx);
    Waveform w(waveform, recordLength, true);

    int startSample;
    double peakAmp;
    bool isNeutron = w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                                startSample, peakAmp);
    if (!isNeutron) {
        std::cerr << "Event " << idx << " did not pass isNeutronLike(); skipping" << std::endl;
        return;
    }

    PSD::Result r = psd.analyzeAt(w.getSubtractedADC(), recordLength, startSample);

    char label[64];
    snprintf(label, sizeof(label), "PSD = %.3f", r.psdRatio);

    overlay.draw(w, startSample, SAMPLE_SPACING_NS, idx, c, pdfName, label);

    std::cout << "Event " << idx << ": startSample=" << startSample
              << ", peakAmp=" << peakAmp << ", PSD=" << r.psdRatio << std::endl;
}

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root> [event1] [event2] ..." << std::endl;
        return 1;
    }

    std::string inFile = Utilities::INPUT_DIR + argv[1];

    std::vector<Long64_t> requestedIndices;
    for (int i = 2; i < argc; i++) requestedIndices.push_back(std::atoll(argv[i]));

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

    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, SAMPLE_SPACING_NS, 0.0, START_SHIFT_NS);
    NeutronCurveOverlay overlay;

    std::string outDir = Utilities::makeOutputDir(argv[0], argv[1]);
    gSystem->mkdir(outDir.c_str(), true);
    std::string pdfName = outDir + "curveOverlay.pdf";

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetGrid();
    c->Print((pdfName + "[").c_str());

    if (!requestedIndices.empty()) {
        // Specific events requested.
        for (Long64_t idx : requestedIndices) {
            if (idx < 0 || idx >= nEntries) {
                std::cerr << "Skipping event " << idx << ": out of range (0-" << nEntries - 1 << ")" << std::endl;
                continue;
            }
            plotOneEvent(tree, idx, recordLength, waveform, psd, overlay, c, pdfName);
        }
    } else {
        // No indices given: scan the whole file automatically.
        int plotted = 0;
        for (Long64_t i = 0; i < nEntries && plotted < MAX_SAMPLE_PLOTS; i++) {
            tree->GetEntry(i);
            Waveform w(waveform, recordLength, true);

            int startSample;
            double peakAmp;
            if (!w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                           startSample, peakAmp)) continue;

            PSD::Result r = psd.analyzeAt(w.getSubtractedADC(), recordLength, startSample);
            char label[64];
            snprintf(label, sizeof(label), "PSD = %.3f", r.psdRatio);

            overlay.draw(w, startSample, SAMPLE_SPACING_NS, i, c, pdfName, label);
            plotted++;
        }
        std::cout << "Scanned file, plotted " << plotted << " neutron-like events (capped at "
                  << MAX_SAMPLE_PLOTS << ")" << std::endl;
    }

    c->Print((pdfName + "]").c_str());
    std::cout << "Saved: " << pdfName << std::endl;

    f->Close();
    return 0;
}