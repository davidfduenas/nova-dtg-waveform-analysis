// selectNeutronBox2D.cpp
// Usage: ./selectNeutronBox2D <file.root>
//
// Applies a 2D (PSD ratio, amplitude) box selection, but ONLY to events
// already classified neutron-like by Waveform::isNeutronLike() -- this
// is a refinement of Sample #1, not a new independent classifier.
//
// Box: PSD in [PSD_MIN, PSD_MAX], amplitude > AMPLITUDE_MIN.
//   PASS   -> selectedNeutron_<file>.pdf
//   REJECT -> rejectedNeutron_<file>.pdf
// Both plotted with gate lines via PSD::drawEventWaveformAt (reusing
// the existing, already-correct detected-pulse anchor).

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

// The 2D box, applied only to neutron-like events.
static const double PSD_MIN        = 0.7;
static const double PSD_MAX        = 1.0;
static const double AMPLITUDE_MIN  = 400.0;

static const int MAX_EVENTS = 1000;

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
    std::string passPdf   = OUTPUT_PATH + "selectedNeutron_" + outName + ".pdf";
    std::string rejectPdf = OUTPUT_PATH + "rejectedNeutron_" + outName + ".pdf";

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetGrid();
    c->Print((passPdf + "[").c_str());
    c->Print((rejectPdf + "[").c_str());

    Long64_t nEntries = tree->GetEntries();
    long nNeutronLike = 0, nPass = 0, nReject = 0;
    int passPlotted = 0, rejectPlotted = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);

        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;

        int startSample;
        double localPeak;
        bool isNeutron = w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                                    startSample, localPeak);
        if (!isNeutron) continue; // only applying this box to neutron-like events

        nNeutronLike++;

        const std::vector<double>& shape = w.getSubtractedADC();
        PSD::Result r = psd.analyzeAt(shape, recordLength, startSample);

        bool inBox = (r.psdRatio >= PSD_MIN && r.psdRatio <= PSD_MAX && localPeak > AMPLITUDE_MIN);

        if (inBox) {
            nPass++;
            if (passPlotted < MAX_EVENTS) {
                psd.drawEventWaveformAt(shape, recordLength, startSample, i, c, passPdf);
                passPlotted++;
            }
        } else {
            nReject++;
            if (rejectPlotted < MAX_EVENTS) {
                psd.drawEventWaveformAt(shape, recordLength, startSample, i, c, rejectPdf);
                rejectPlotted++;
            }
        }
    }

    c->Print((passPdf + "]").c_str());
    c->Print((rejectPdf + "]").c_str());

    std::cout << "Neutron-like events: " << nNeutronLike << std::endl;
    std::cout << "Pass 2D box (PSD [" << PSD_MIN << "," << PSD_MAX << "], amplitude > " << AMPLITUDE_MIN
              << "): " << nPass << " (" << 100.0 * nPass / nNeutronLike << "%)" << std::endl;
    std::cout << "Rejected: " << nReject << " (" << 100.0 * nReject / nNeutronLike << "%)" << std::endl;
    std::cout << "Plotted " << passPlotted << " passed -> " << passPdf << std::endl;
    std::cout << "Plotted " << rejectPlotted << " rejected -> " << rejectPdf << std::endl;

    f->Close();
    return 0;
}