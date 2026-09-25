// inspectPSDRange.cpp
// Usage: ./inspectPSDRange <file.root>
//
// Dumps waveforms + all PSD values for events whose PSD ratio falls in
// [PSD_MIN, PSD_MAX]. Meant for looking at the suspicious high-PSD
// events directly.

#include "Waveform.h"
#include "PSD.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TError.h>

#include <iostream>
#include <string>

static const std::string OUTPUT_PATH = "results/";
static const std::string INPUT_DIR   = "/Users/david/DTGAnalysis/data/testruns/";

static const double POST_TRIGGER_PERCENT = 80.0;
static const double N_COEFFICIENT        = 8.0;
static const double CONSTANT_LATENCY     = 73.2889;
static const double NS_PER_SAMPLE        = 2.0;

static const double SHORT_GATE_NS  = 40.0;
static const double LONG_GATE_NS   = 300.0;
static const double START_SHIFT_NS = 10.0;

static const double AMPLITUDE_CUTOFF = 50.0;

// PSD range to inspect
static const double PSD_MIN = 0.4;
static const double PSD_MAX = 1.0;

static const int MAX_EVENTS = 1300;

// Control window: sums the same number of samples as the long gate, but
// from far away from the pulse (pure baseline, no signal). If this comes
// out close to (qLong - qShort), that confirms the excess in qLong is a
// baseline bias present throughout the record, not something specific to
// the region right after the pulse.
static const int CONTROL_START_SAMPLE = 600; // ~1200 ns -- far from any pulse

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning; // suppress ROOT's "Info in <TCanvas::Print>" spam

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
            SHORT_GATE_NS, LONG_GATE_NS, NS_PER_SAMPLE, AMPLITUDE_CUTOFF,
            START_SHIFT_NS);

    std::string outName = argv[1];
    for (char& ch : outName) if (ch == '.') ch = '_';
    std::string pdfName = OUTPUT_PATH + "inspectPSD_" + outName + ".pdf";

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetGrid();
    gErrorIgnoreLevel = kWarning; // set again here -- ROOT can reset this during init
    c->Print((pdfName + "[").c_str());

    std::cout << "PSD range: [" << PSD_MIN << ", " << PSD_MAX << "]"
              << "   amplitude cutoff: " << AMPLITUDE_CUTOFF << std::endl;
    std::cout << "\nevent\tmaxAmp\tqShort\t\tqLong\t\tPSD\t\tstart\tshortEnd\tlongEnd\t\tcontrolSum" << std::endl;

    Long64_t nEntries = tree->GetEntries();
    int found = 0;

    for (Long64_t i = 0; i < nEntries && found < MAX_EVENTS; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);

        if (w.getMaxAmpADC() <= AMPLITUDE_CUTOFF) continue;

        PSD::Result r = psd.analyze(w.getSubtractedADC(), recordLength);
        if (r.psdRatio < PSD_MIN || r.psdRatio > PSD_MAX) continue;

        // Control sum: same width as the long gate (longEnd - start
        // samples), taken from a quiet region far from the pulse.
        const std::vector<double>& shape = w.getSubtractedADC();
        int controlWidth = r.longGateEndSample - r.pulseStartSample;
        int controlEnd = std::min(CONTROL_START_SAMPLE + controlWidth, recordLength);
        double controlSum = 0;
        for (int s = CONTROL_START_SAMPLE; s < controlEnd; s++) controlSum += shape[s];

        std::cout << i << "\t" << w.getMaxAmpADC()
                  << "\t" << r.qShort
                  << "\t" << r.qLong
                  << "\t" << r.psdRatio
                  << "\t" << r.pulseStartSample
                  << "\t" << r.shortGateEndSample
                  << "\t\t" << r.longGateEndSample
                  << "\t\t" << controlSum
                  << "\t(qLong-qShort=" << (r.qLong - r.qShort) << ")" << std::endl;

        psd.drawEventWaveform(w.getSubtractedADC(), recordLength, i, c, pdfName);
        found++;
    }

    c->Print((pdfName + "]").c_str());

    std::cout << "\nFound " << found << " events in range (capped at " << MAX_EVENTS << ")" << std::endl;
    std::cout << "Saved: " << pdfName << std::endl;

    f->Close();
    return 0;
}