// checkDecliningEnvelope.cpp
// Usage: ./checkDecliningEnvelope <file.root>
// Checks two known events against hasDecliningEnvelope(), then scans the
// full neutron-like population and plots declining vs flat events.

#include "Waveform.h"
#include "PSD.h"
#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <cstdio>

static const double SAMPLE_SPACING_NS = 2.0;
static const double AMPLITUDE_CUTOFF = 50.0;
static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;

static const int NUM_CHUNKS = 5;
static const int MIN_DECLINING_PAIRS = 3;

// Two known events to check specifically -- update to match the file under test.
static const long KNOWN_EVENT_1 = 653149;
static const long KNOWN_EVENT_2 = 107670;

static const int MAX_SAMPLE_PLOTS = 200;

static const double POST_TRIGGER_PERCENT = 80.0;
static const double N_COEFFICIENT        = 8.0;
static const double CONSTANT_LATENCY     = 73.2889;
static const double SHORT_GATE_NS        = 40.0;
static const double LONG_GATE_NS         = 300.0;
static const double START_SHIFT_NS       = 7.0;

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input_file.root>" << std::endl;
        return 1;
    }

    std::string inFile = Utilities::INPUT_DIR + argv[1];

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
    std::cout << inFile << ": " << nEntries << " total events" << std::endl;

    std::cout << "\n--- Specific known events ---" << std::endl;
    for (long idx : {KNOWN_EVENT_1, KNOWN_EVENT_2}) {
        if (idx < 0 || idx >= nEntries) {
            std::cout << "Event " << idx << ": out of range for this file" << std::endl;
            continue;
        }
        tree->GetEntry(idx);
        Waveform w(waveform, recordLength, true);

        int startSample; double peakAmp;
        bool isNeutron = w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                                    startSample, peakAmp);
        bool declining = isNeutron && w.hasDecliningEnvelope(startSample, NUM_CHUNKS, MIN_DECLINING_PAIRS);

        std::cout << "Event " << idx << ": isNeutronLike=" << isNeutron
                  << "  hasDecliningEnvelope=" << declining << std::endl;
    }

    std::cout << "\n--- Full neutron-like population ---" << std::endl;

    std::string outDir = Utilities::makeOutputDir(argv[0], argv[1]);
    gSystem->mkdir(outDir.c_str(), true);
    std::string decliningPdf = outDir + "declining.pdf";
    std::string flatPdf      = outDir + "flat.pdf";

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetGrid();
    c->Print((decliningPdf + "[").c_str());
    c->Print((flatPdf + "[").c_str());

    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, SAMPLE_SPACING_NS, 0.0, START_SHIFT_NS);

    long nNeutronLike = 0, nDeclining = 0, nFlat = 0;
    int decliningPlotted = 0, flatPlotted = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);

        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;

        int startSample; double peakAmp;
        bool isNeutron = w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                                    startSample, peakAmp);
        if (!isNeutron) continue;

        nNeutronLike++;
        bool declining = w.hasDecliningEnvelope(startSample, NUM_CHUNKS, MIN_DECLINING_PAIRS);

        if (declining) {
            nDeclining++;
            if (decliningPlotted < MAX_SAMPLE_PLOTS) {
                const std::vector<double>& shape = w.getSubtractedADC();
                PSD::Result r = psd.analyzeAt(shape, recordLength, startSample);

                char label[64];
                snprintf(label, sizeof(label), "#splitline{PSD = %.3f}{Declining}", r.psdRatio);

                w.drawWithThreshold(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES, i,
                                     SAMPLE_SPACING_NS, c, decliningPdf, label);
                decliningPlotted++;
            }
        } else {
            nFlat++;
            if (flatPlotted < MAX_SAMPLE_PLOTS) {
                const std::vector<double>& shape = w.getSubtractedADC();
                PSD::Result r = psd.analyzeAt(shape, recordLength, startSample);

                char label[64];
                snprintf(label, sizeof(label), "#splitline{PSD = %.3f}{Flat / no decline}", r.psdRatio);

                w.drawWithThreshold(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES, i,
                                     SAMPLE_SPACING_NS, c, flatPdf, label);
                flatPlotted++;
            }
        }
    }

    c->Print((decliningPdf + "]").c_str());
    c->Print((flatPdf + "]").c_str());

    std::cout << "Total neutron-like: " << nNeutronLike << std::endl;
    std::cout << "Declining (likely real): " << nDeclining << " (" << 100.0 * nDeclining / nNeutronLike << "%)" << std::endl;
    std::cout << "Flat/no decline (likely false positive): " << nFlat << " (" << 100.0 * nFlat / nNeutronLike << "%)" << std::endl;
    std::cout << "Plotted " << decliningPlotted << " declining -> " << decliningPdf << std::endl;
    std::cout << "Plotted " << flatPlotted << " flat -> " << flatPdf << std::endl;

    f->Close();
    return 0;
}
