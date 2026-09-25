// countNeutronLikeByAmplitude.cpp
// Usage: ./countNeutronLikeByAmplitude <file.root>
//
// Finds every event classified neutron-like by Waveform::isNeutronLike(),
// then prints how many of them have amplitude above each threshold in
// AMPLITUDE_THRESHOLDS -- shows directly how the total neutron-like count
// changes as the minimum amplitude requirement is raised.

#include "Waveform.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>

static const std::string INPUT_DIR = "/Users/david/DTGAnalysis/data/testruns/";

static const double AMPLITUDE_CUTOFF = 50.0; // minimum to even be a "candidate" at all
static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;

// Minimum-amplitude thresholds to scan -- edit this list as needed.
static const std::vector<double> AMPLITUDE_THRESHOLDS = {
    0, 50, 100, 200, 300, 400, 500, 600, 800, 1000, 1500, 2000, 3000, 4000, 5000
};

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

    // Collect the local peak amplitude of every neutron-like event once.
    std::vector<double> neutronLikeAmps;

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);

        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;

        int startSample;
        double localPeak;
        bool isNeutron = w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD, NEUTRON_MIN_STABLE_SAMPLES,
                                                    startSample, localPeak);
        if (!isNeutron) continue;

        neutronLikeAmps.push_back(localPeak);
    }

    std::cout << "Total neutron-like events (no minimum beyond " << AMPLITUDE_CUTOFF
              << " ADC candidate cutoff): " << neutronLikeAmps.size() << std::endl;

    std::cout << "\nMinimum amplitude (ADC)\tCount above it" << std::endl;
    for (double thr : AMPLITUDE_THRESHOLDS) {
        long count = 0;
        for (double amp : neutronLikeAmps) if (amp > thr) count++;
        std::cout << thr << "\t\t\t" << count << std::endl;
    }

    f->Close();
    return 0;
}