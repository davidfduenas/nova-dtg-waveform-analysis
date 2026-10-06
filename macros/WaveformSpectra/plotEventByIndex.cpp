// plotEventByIndex.cpp
// Usage: ./plotEventByIndex <file.root> <event1> [event2] [event3] ...
// Plots the full waveform (baseline-subtracted and raw) for specific
// event indices, no amplitude filtering.

#include "Waveform.h"
#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

static const double SAMPLE_SPACING_NS = 2.0;

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc < 3) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_file.root> <event1> [event2] [event3] ..." << std::endl;
        return 1;
    }

    std::string inFile = Utilities::INPUT_DIR + argv[1];
    std::vector<Long64_t> eventIndices;
    for (int i = 2; i < argc; i++) {
        eventIndices.push_back(std::atoll(argv[i]));
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
    tree->SetBranchAddress("recordLength", &recordLength);
    tree->SetBranchAddress("waveform", waveform);

    Long64_t nEntries = tree->GetEntries();
    std::cout << inFile << ": " << nEntries << " total events" << std::endl;

    std::string outDir = Utilities::makeOutputDir(argv[0], argv[1]);
    gSystem->mkdir(outDir.c_str(), true);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    TCanvas* c = new TCanvas("c", "Event Waveforms", 900, 700);
    c->SetGrid();

    std::string pdfName = outDir + "events_by_index.pdf";
    c->Print((pdfName + "[").c_str());

    for (Long64_t idx : eventIndices) {
        if (idx < 0 || idx >= nEntries) {
            std::cerr << "Skipping event " << idx << ": out of range (0-"
                      << nEntries - 1 << ")" << std::endl;
            continue;
        }

        tree->GetEntry(idx);
        Waveform w(waveform, recordLength, true);

        const std::vector<double>& shape = w.getSubtractedADC();

        TH1D* h = new TH1D(("h" + std::to_string(idx)).c_str(),
                            ("Event " + std::to_string(idx)).c_str(),
                            recordLength, 0, recordLength * SAMPLE_SPACING_NS);
        for (int s = 0; s < recordLength; s++) {
            h->SetBinContent(s + 1, shape[s]);
        }

        h->GetXaxis()->SetTitle("Time (ns)");
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->SetTitle("Baseline - Signal (ADC)");
        h->GetYaxis()->CenterTitle();
        h->SetLineColor(kBlue);

        h->Draw("HIST");
        c->Print(pdfName.c_str());

        std::cout << "Event " << idx << ": maxAmp=" << w.getMaxAmpADC() << " ADC" << std::endl;

        // Raw (non-baseline-subtracted) waveform, to check whether the
        // pretrigger region is genuinely flat before baseline calculation.
        TH1D* hRaw = new TH1D(("hRaw" + std::to_string(idx)).c_str(),
                               ("Event " + std::to_string(idx) + " (raw)").c_str(),
                               recordLength, 0, recordLength * SAMPLE_SPACING_NS);
        for (int s = 0; s < recordLength; s++) {
            hRaw->SetBinContent(s + 1, waveform[s]);
        }
        hRaw->GetXaxis()->SetTitle("Time (ns)");
        hRaw->GetXaxis()->CenterTitle();
        hRaw->GetYaxis()->SetTitle("Raw ADC");
        hRaw->GetYaxis()->CenterTitle();
        hRaw->SetLineColor(kBlack);

        hRaw->Draw("HIST");
        c->Print(pdfName.c_str());

        int nCheck = std::min(100, recordLength);
        double sum = 0, sum2 = 0;
        for (int s = 0; s < nCheck; s++) { sum += waveform[s]; sum2 += (double)waveform[s]*waveform[s]; }
        double mean = sum / nCheck;
        double var = (sum2 - nCheck*mean*mean) / (nCheck - 1);
        double rms = (var > 0) ? std::sqrt(var) : 0;
        std::cout << "  First " << nCheck << " raw samples: mean=" << mean
                  << ", rms=" << rms << std::endl;

        delete h;
        delete hRaw;
    }

    c->Print((pdfName + "]").c_str());
    std::cout << "Saved: " << pdfName << std::endl;

    f->Close();
    return 0;
}
