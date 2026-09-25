// showGateLocations.cpp
// Usage: ./showGateLocations <file.root>
//
// Classifies each candidate event the same way as runPSDSamples13.cpp,
// but anchors PSD's gates at the ACTUAL detected pulse location (from
// Waveform::findNeutronLikeStretch() / countSeparatePulses()) instead of
// the fixed acquisition-timing formula. Amplitude used for PSD's x-axis
// is also the LOCAL peak of that same detected pulse, not the whole-
// record max -- so both axes always come from the same physical event.
//
// Produces two annotated PDFs (gate lines drawn via PSD::drawEventWaveformAt):
//   gateLocations_neutronLike_<file>.pdf
//   gateLocations_gammaLike_<file>.pdf

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

static const double SHORT_GATE_NS = 40.0;
static const double LONG_GATE_NS  = 300.0;
static const double START_SHIFT_NS = 7.0; // pre-gate, applied uniformly to both samples

static const double AMPLITUDE_CUTOFF = 50.0;

// isNeutronLike() / countSeparatePulses() parameters.
static const double NEUTRON_AMPLITUDE_THRESHOLD = 50.0;
static const int    NEUTRON_MIN_STABLE_SAMPLES  = 38;
static const int    PILEUP_MIN_QUIET_SAMPLES    = 15;

static const int MAX_SAMPLE_PLOTS = 200;

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

    // Note: gate widths are fixed regardless of anchor -- only the START
    // point changes (explicit detected location instead of the formula).
    PSD psd(POST_TRIGGER_PERCENT, N_COEFFICIENT, CONSTANT_LATENCY,
            SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, 0.0, START_SHIFT_NS);

    std::string outName = argv[1];
    for (char& ch : outName) if (ch == '.') ch = '_';
    std::string neutronPdf = OUTPUT_PATH + "gateLocations_neutronLike_" + outName + ".pdf";
    std::string gammaPdf   = OUTPUT_PATH + "gateLocations_gammaLike_" + outName + ".pdf";

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetGrid();
    c->Print((neutronPdf + "[").c_str());
    c->Print((gammaPdf + "[").c_str());

    Long64_t nEntries = tree->GetEntries();
    long nCandidates = 0, nNeutronLike = 0, nGammaLike = 0, nPileup = 0;
    int neutronPlotted = 0, gammaPlotted = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);

        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;
        nCandidates++;

        const std::vector<double>& shape = w.getSubtractedADC();

        int startSample;
        double localPeak;
        bool isNeutron = w.findNeutronLikeStretch(NEUTRON_AMPLITUDE_THRESHOLD,
                                                    NEUTRON_MIN_STABLE_SAMPLES,
                                                    startSample, localPeak);

        if (isNeutron) {
            nNeutronLike++;
            if (neutronPlotted < MAX_SAMPLE_PLOTS) {
                psd.drawEventWaveformAt(shape, recordLength, startSample, i, c, neutronPdf);
                neutronPlotted++;
            }
            continue;
        }

        std::vector<int> pulseStarts;
        auto pulses = w.countSeparatePulses(NEUTRON_AMPLITUDE_THRESHOLD, PILEUP_MIN_QUIET_SAMPLES, pulseStarts);

        if (pulses.size() == 1) {
            nGammaLike++;
            if (gammaPlotted < MAX_SAMPLE_PLOTS) {
                psd.drawEventWaveformAt(shape, recordLength, pulseStarts[0], i, c, gammaPdf);
                gammaPlotted++;
            }
        } else {
            nPileup++;
        }
    }

    c->Print((neutronPdf + "]").c_str());
    c->Print((gammaPdf + "]").c_str());

    std::cout << "Total candidates: " << nCandidates << std::endl;
    std::cout << "Neutron-like: " << nNeutronLike << std::endl;
    std::cout << "Gamma-like (single pulse): " << nGammaLike << std::endl;
    std::cout << "Pile-up (excluded): " << nPileup << std::endl;
    std::cout << "Plotted " << neutronPlotted << " neutron-like -> " << neutronPdf << std::endl;
    std::cout << "Plotted " << gammaPlotted << " gamma-like -> " << gammaPdf << std::endl;

    f->Close();
    return 0;
}