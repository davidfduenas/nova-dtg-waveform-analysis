// plotWaveforms.cpp
// Scans events in a converted ROOT file, finds the first N events whose
// max amplitude falls within [low, high] (in ADC counts or mV, selectable),
// and saves each one's baseline-subtracted waveform as a page in one PDF.
// The PDF's first page shows the run duration and event count. The output
// filename includes the record length and the selected amplitude range.
//
// Usage:
//   ./plotWaveforms <input_file.root> <low> <high> <N> [unit: adc|mv]

#include "Waveform.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPaveText.h>

#include <iostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

std::string detectLabel(const std::string& path)
{
    std::string lower = path;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    if (lower.find("bkg") != std::string::npos) return "bkg";
    return "dtgon";
}

// Formats a number for use in a filename: whole numbers with no decimal
// point, non-whole numbers with one decimal place -- avoids ugly trailing
// zeros like "200.00" while still handling fractional mV values cleanly.
std::string formatAmpForFilename(double val)
{
    std::ostringstream oss;
    if (std::fabs(val - std::round(val)) < 1e-6) {
        oss << (long long)std::round(val);
    } else {
        oss.precision(1);
        oss << std::fixed << val;
    }
    return oss.str();
}

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    static const std::string DATA_PATH   = "../data/testruns/";
    static const std::string OUTPUT_PATH = "results/";
    static const double SAMPLE_SPACING_NS = 2.0;

    if (argc < 5 || argc > 6) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_file.root> <low> <high> <N> [unit: adc|mv]" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    double low  = std::atof(argv[2]);
    double high = std::atof(argv[3]);
    int N = std::atoi(argv[4]);
    std::string unit = (argc == 6) ? argv[5] : "adc";

    std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);
    if (unit != "adc" && unit != "mv") {
        std::cerr << "ERROR: unit must be 'adc' or 'mv'" << std::endl;
        return 1;
    }
    bool useMv = (unit == "mv");
    std::string unitLabel = useMv ? "mV" : "ADC";

    std::string label = detectLabel(inputFile);

    std::string fullPath = DATA_PATH + inputFile;
    TFile* f = TFile::Open(fullPath.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open " << fullPath << std::endl;
        return 1;
    }

    TTree* tree = (TTree*)f->Get("waveforms");
    if (!tree) {
        std::cerr << "ERROR: could not find tree 'waveforms' in " << fullPath << std::endl;
        return 1;
    }

    Int_t recordLength;
    int waveform[20000];
    Double_t unixTime;
    tree->SetBranchAddress("recordLength", &recordLength);
    tree->SetBranchAddress("waveform", waveform);
    tree->SetBranchAddress("unixTime", &unixTime);

    Long64_t nEntries = tree->GetEntries();

    tree->GetEntry(0);
    double t_start = unixTime;
    tree->GetEntry(nEntries - 1);
    double t_end = unixTime;
    double duration = t_end - t_start;

    gSystem->mkdir(OUTPUT_PATH.c_str(), true);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    int recNs = std::round(recordLength * SAMPLE_SPACING_NS);
    int recUs = std::round(recNs / 1000.0);

    std::string rangeLabel = formatAmpForFilename(low) + "to" + formatAmpForFilename(high) + unitLabel;
    std::string pdfName = OUTPUT_PATH + "waveforms_" + label + "_" + std::to_string(recUs) + "us_"
                         + rangeLabel + ".pdf";

    TCanvas* c = new TCanvas("c", "", 900, 700);
    c->SetGrid();

    c->Print((pdfName + "[").c_str());

    // --- Info page ---
    c->Clear();
    c->SetGrid(0,0);
    TPaveText* info = new TPaveText(0.1, 0.3, 0.9, 0.8, "NDC");
    info->SetBorderSize(0);
    info->SetFillStyle(0);
    info->SetTextAlign(12);
    info->AddText(("File: " + inputFile).c_str());
    info->AddText(Form("Total events: %lld", nEntries));
    info->AddText(Form("Run duration: %.3f s", duration));
    info->AddText(Form("Record length: %d samples (%d us)", (int)recordLength, recUs));
    info->AddText(Form("Selection range: %.2f - %.2f %s", low, high, unitLabel.c_str()));
    info->Draw();
    c->Print(pdfName.c_str());
    c->SetGrid();

    int found = 0;
    for (Long64_t i = 0; i < nEntries && found < N; i++) {
        tree->GetEntry(i);

        Waveform w(waveform, recordLength, true);
        double amp = useMv ? w.getMaxAmpVolts() : w.getMaxAmpADC();

        if (amp < low || amp > high) continue;

        found++;

        std::vector<double> shape = useMv ? w.getSubtractedVolts() : w.getSubtractedADC();
        std::string title = "Event " + std::to_string(i);
        TH1D* h = new TH1D(Form("wf_%d", found), title.c_str(), recordLength,
                            0, recordLength * SAMPLE_SPACING_NS);
        for (int s = 0; s < recordLength; s++) h->SetBinContent(s+1, shape[s]);

        h->SetLineColor(kBlue);
        h->SetLineWidth(2);
        h->GetXaxis()->SetTitle("Time (ns)");
        h->GetXaxis()->CenterTitle();
        h->GetYaxis()->SetTitle(("Baseline - Signal (" + unitLabel + ")").c_str());
        h->GetYaxis()->CenterTitle();
        h->Draw("HIST");

        c->Print(pdfName.c_str());

        std::cout << "Page " << found+1 << ": event " << i
                  << ", maxAmp=" << amp << " " << unitLabel << std::endl;

        delete h;
    }

    c->Print((pdfName + "]").c_str());

    if (found < N) {
        std::cout << "\nOnly found " << found << " events in range [" << low << ", " << high
                  << "] " << unitLabel << " (requested " << N << ")" << std::endl;
    }

    std::cout << "\nSaved: " << pdfName << " (" << found+1 << " pages, including info page)" << std::endl;

    f->Close();
    return 0;
}