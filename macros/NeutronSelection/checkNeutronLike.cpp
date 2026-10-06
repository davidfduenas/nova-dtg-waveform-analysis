// checkNeutronLike.cpp
// Usage: ./checkNeutronLike <file.root>
// Classifies each candidate event as neutron-like, gamma-like (single
// pulse), or pile-up, and plots a sample of each.

#include "Waveform.h"
#include "PSD.h"
#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>

static const double SAMPLE_SPACING_NS = 2.0;
static const double AMPLITUDE_CUTOFF = 50.0;

static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;
static const int    PILEUP_MIN_QUIET_SAMPLES    = 15;

static const int MAX_SAMPLE_PLOTS = 1000;

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

    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    std::string outDir = Utilities::makeOutputDir(argv[0], argv[1]);
    gSystem->mkdir(outDir.c_str(), true);

    std::string neutronLikePdf = outDir + "neutronLike.pdf";
    std::string gammaLikePdf   = outDir + "gammaLike.pdf";
    std::string pileupPdf      = outDir + "pileup.pdf";

    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, SAMPLE_SPACING_NS, 0.0, START_SHIFT_NS);

    long nCandidates = 0, nNeutronLike = 0, nGammaLike = 0, nPileup = 0;
    int neutronLikePlotted = 0, gammaLikePlotted = 0, pileupPlotted = 0;

    TCanvas* c2 = new TCanvas("c2", "", 900, 700);
    c2->SetGrid();
    c2->Print((neutronLikePdf + "[").c_str());
    c2->Print((gammaLikePdf + "[").c_str());
    c2->Print((pileupPdf + "[").c_str());

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);
        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;
        nCandidates++;

        const std::vector<double>& shape = w.getSubtractedADC();

        bool result = w.isNeutronLike(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES);
        if (result) {
            nNeutronLike++;
            if (neutronLikePlotted < MAX_SAMPLE_PLOTS) {
                int startSample; double peakAmp;
                w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                          startSample, peakAmp);
                PSD::Result r = psd.analyzeAt(shape, recordLength, startSample);

                char label[64];
                snprintf(label, sizeof(label), "PSD = %.3f", r.psdRatio);

                w.drawWithThreshold(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES, i,
                                     SAMPLE_SPACING_NS, c2, neutronLikePdf, label);
                neutronLikePlotted++;
            }
        } else {
            std::vector<int> pulseStarts;
            auto pulses = w.countSeparatePulses(NEUTRON_AMPLITUDE_THRESHOLD, PILEUP_MIN_QUIET_SAMPLES, pulseStarts);
            if (pulses.size() == 1) {
                nGammaLike++;
                if (gammaLikePlotted < MAX_SAMPLE_PLOTS) {
                    PSD::Result r = psd.analyzeAt(shape, recordLength, pulseStarts[0]);

                    char label[64];
                    snprintf(label, sizeof(label), "PSD = %.3f", r.psdRatio);

                    w.drawWithThreshold(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES, i,
                                         SAMPLE_SPACING_NS, c2, gammaLikePdf, label);
                    gammaLikePlotted++;
                }
            } else {
                nPileup++;
                if (pileupPlotted < MAX_SAMPLE_PLOTS) {
                    int anchorStart = pulseStarts.empty() ? 0 : pulseStarts[0];
                    PSD::Result r = psd.analyzeAt(shape, recordLength, anchorStart);

                    char label[64];
                    snprintf(label, sizeof(label), "PSD = %.3f", r.psdRatio);

                    w.drawWithThreshold(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES, i,
                                         SAMPLE_SPACING_NS, c2, pileupPdf, label);
                    pileupPlotted++;
                }
            }
        }
    }

    c2->Print((neutronLikePdf + "]").c_str());
    c2->Print((gammaLikePdf + "]").c_str());
    c2->Print((pileupPdf + "]").c_str());

    std::cout << "Total candidates: " << nCandidates << std::endl;
    std::cout << "Neutron-like: " << nNeutronLike << " / " << nCandidates
              << " (" << 100.0 * nNeutronLike / nCandidates << "%)" << std::endl;
    std::cout << "Gamma-pulse-like (single pulse): " << nGammaLike << " / " << nCandidates
              << " (" << 100.0 * nGammaLike / nCandidates << "%)" << std::endl;
    std::cout << "Pile-up (multi pulse): " << nPileup << " / " << nCandidates
              << " (" << 100.0 * nPileup / nCandidates << "%)" << std::endl;
    std::cout << "Plotted " << neutronLikePlotted << " neutron-like samples -> " << neutronLikePdf << std::endl;
    std::cout << "Plotted " << gammaLikePlotted << " gamma-pulse-like samples -> " << gammaLikePdf << std::endl;
    std::cout << "Plotted " << pileupPlotted << " pile-up samples -> " << pileupPdf << std::endl;

    f->Close();
    return 0;
}
