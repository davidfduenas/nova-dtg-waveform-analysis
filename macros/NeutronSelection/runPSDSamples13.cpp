// runPSDSamples13.cpp
// Usage: ./runPSDSamples13 <file.root>
//
// Classifies every candidate event (amplitude > AMPLITUDE_CUTOFF) into:
//   Sample #1: isNeutronLike() == true                 -> neutron candidates
//   Sample #2: isNeutronLike() == false                 -> gammas + pile-up
//     -> split further via countSeparatePulses():
//        single pulse (size==1)  -> Sample #3 (clean single gammas)
//        multi pulse  (size>=2)  -> excluded (pile-up)
//
// Runs PSD on Sample #1 and Sample #3 separately, saving two independent
// sets of plots.
//
// IMPORTANT: POST_TRIGGER_PERCENT/CONSTANT_LATENCY/gate widths below must
// match how this file was actually acquired -- see PSD.h for details.

#include "Waveform.h"
#include "PSD.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>

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

// isNeutronLike() / countSeparatePulses() parameters.
static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;
static const int    PILEUP_MIN_QUIET_SAMPLES    = 15;

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

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

    // amplitudeCutoff=0 for both PSD instances: selection into each
    // sample already happened via isNeutronLike()/countSeparatePulses()
    // before addEvent() is called, so PSD shouldn't re-filter by amplitude.
    PSD psdSample1(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
                   SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, 0.0, START_SHIFT_NS);
    PSD psdSample3(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
                   SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, 0.0, START_SHIFT_NS);

    Long64_t nEntries = tree->GetEntries();
    std::cout << inFile << ": " << nEntries << " total events" << std::endl;

    long nCandidates = 0;
    long nSample1 = 0, nSample2 = 0, nSample3 = 0, nPileupExcluded = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);

        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;
        nCandidates++;

        const std::vector<double>& shape = w.getSubtractedADC();

        int startSample;
        double localPeak;
        if (w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                      startSample, localPeak)) {
            nSample1++;
            psdSample1.addEventAt(shape, recordLength, localPeak, startSample);
        } else {
            nSample2++;
            std::vector<int> pulseStarts;
            auto pulses = w.countSeparatePulses(NEUTRON_AMPLITUDE_THRESHOLD, PILEUP_MIN_QUIET_SAMPLES, pulseStarts);
            if (pulses.size() == 1) {
                nSample3++;
                psdSample3.addEventAt(shape, recordLength, pulses[0], pulseStarts[0]);
            } else {
                nPileupExcluded++;
            }
        }
    }

    std::string outName = argv[1];
    for (char& ch : outName) if (ch == '.') ch = '_';

    psdSample1.savePlots(OUTPUT_PATH, "sample1_neutronlike_" + outName);
    psdSample3.savePlots(OUTPUT_PATH, "sample3_singlegamma_" + outName);
    PSD::saveCombinedPlot(psdSample3, psdSample1, OUTPUT_PATH, outName);

    std::cout << "\nTotal candidates (amplitude > " << AMPLITUDE_CUTOFF << "): " << nCandidates << std::endl;
    std::cout << "Sample #1 (neutron-like): " << nSample1
               << " (" << 100.0 * nSample1 / nCandidates << "%)" << std::endl;
    std::cout << "Sample #2 (not neutron-like): " << nSample2
               << " (" << 100.0 * nSample2 / nCandidates << "%)" << std::endl;
    std::cout << "  -> Sample #3 (single pulse, kept): " << nSample3
               << " (" << 100.0 * nSample3 / nCandidates << "% of total)" << std::endl;
    std::cout << "  -> Pile-up (multi pulse, excluded): " << nPileupExcluded
               << " (" << 100.0 * nPileupExcluded / nCandidates << "% of total)" << std::endl;

    std::cout << "\nSaved PSD plots for Sample #1 -> psdRatio_sample1_neutronlike_" << outName << ".pdf"
               << " and psdVsAmplitude_sample1_neutronlike_" << outName << ".pdf" << std::endl;
    std::cout << "Saved PSD plots for Sample #3 -> psdRatio_sample3_singlegamma_" << outName << ".pdf"
               << " and psdVsAmplitude_sample3_singlegamma_" << outName << ".pdf" << std::endl;
    std::cout << "Saved combined plot -> psdVsAmplitude_combined_" << outName << ".pdf" << std::endl;

    f->Close();
    return 0;
}