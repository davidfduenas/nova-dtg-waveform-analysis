// gammaSpectrum.cpp
// Pulse-height spectrum for a gamma run, with a background run overlaid on
// the same plot. Raw counts only -- no rate/duration normalization, no
// subtraction.
//
// Input directory is hardcoded (see INPUT_DIR below); pass just the
// filenames, not full paths.
//
// Usage:
//   ./gammaSpectrum <gamma_file.root> <bkg_file.root> [unit: adc|mv]

#include "Waveform.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLegend.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

static const int    SPECTRUM_N_BINS = 300;
static const double AMP_MIN_ADC     = 0.0;
static const double AMP_MAX_ADC     = 5000.0;
static const bool   USE_LOG_Y       = true;
static const std::string OUTPUT_PATH = "results/";

// Hardcoded input directory -- pass just filenames on the command line.
static const std::string INPUT_DIR = "/Users/david/DTGAnalysis/data/gammaruns/";

std::string recLabel(int recordLength)
{
    double us = recordLength * 2.0 / 1000.0; // 2 ns/sample
    char buf[32];
    snprintf(buf, sizeof(buf), "%gus", us);
    return std::string(buf);
}

// Reads all events' max amplitude from one file. Returns record length via
// recLenOut. Amplitudes are returned in raw ADC counts (unscaled).
std::vector<double> readAmps(const std::string& filename, int& recLenOut)
{
    TFile* f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open file: " << filename << std::endl;
        std::exit(1);
    }

    TTree* tree = (TTree*)f->Get("waveforms");
    if (!tree) {
        std::cerr << "Could not find tree 'waveforms' in file: " << filename << std::endl;
        std::exit(1);
    }

    Int_t recordLength;
    int waveform[20000];

    tree->SetBranchAddress("recordLength", &recordLength);
    tree->SetBranchAddress("waveform", waveform);

    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "No entries in file: " << filename << std::endl;
        std::exit(1);
    }

    tree->GetEntry(0);
    int recLen = recordLength;

    std::vector<double> amps(nEntries);
    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recLen, false);
        amps[i] = w.getMaxAmpADC();
    }

    f->Close();

    std::cout << filename << ": " << nEntries << " events, recordLength = "
              << recLen << std::endl;

    recLenOut = recLen;
    return amps;
}

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <gamma_file.root> <bkg_file.root> [unit: adc|mv]" << std::endl;
        return 1;
    }

    std::string gammaFile = INPUT_DIR + argv[1];
    std::string bkgFile   = INPUT_DIR + argv[2];
    std::string unit = (argc == 4) ? argv[3] : "adc";

    std::string unitLabel = (unit == "mv") ? "mV" : "ADC";
    double unitScale = (unit == "mv") ? Waveform::MV_PER_COUNT : 1.0;

    int gammaRecLen = 0, bkgRecLen = 0;
    std::vector<double> gammaAmps = readAmps(gammaFile, gammaRecLen);
    std::vector<double> bkgAmps   = readAmps(bkgFile, bkgRecLen);

    double ampMin = AMP_MIN_ADC * unitScale;
    double ampMax = AMP_MAX_ADC * unitScale;

    TH1D* hGamma = new TH1D("hGamma", "Gamma", SPECTRUM_N_BINS, ampMin, ampMax);
    TH1D* hBkg   = new TH1D("hBkg",   "Background", SPECTRUM_N_BINS, ampMin, ampMax);

    for (double a : gammaAmps) hGamma->Fill(a * unitScale);
    for (double a : bkgAmps)   hBkg->Fill(a * unitScale);

    // Scale background to match gamma's total area (equal integrals),
    // expressed in raw counts (not unit area / not rate).
    if (hBkg->Integral() > 0) {
        double localCountRatio = hGamma->Integral() / hBkg->Integral();
        hBkg->Scale(localCountRatio);
    }

    gSystem->mkdir(OUTPUT_PATH.c_str(), true);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    TCanvas* c1 = new TCanvas("c1", "Pulse Height Spectrum", 900, 700);
    c1->SetGrid();
    if (USE_LOG_Y) c1->SetLogy();

    hBkg->SetLineColor(kBlack);
    hBkg->SetLineWidth(2);
    hBkg->GetXaxis()->SetTitle(("Max Amplitude (" + unitLabel + ")").c_str());
    hBkg->GetXaxis()->CenterTitle();
    hBkg->GetYaxis()->SetTitle("Events");
    hBkg->GetYaxis()->CenterTitle();

    hGamma->SetLineColor(kRed);
    hGamma->SetLineWidth(2);

    double maxVal = std::max(hBkg->GetMaximum(), hGamma->GetMaximum());
    hBkg->SetMaximum(maxVal * 1.5);
    hBkg->SetMinimum(0.5);

    hBkg->Draw("HIST");
    hGamma->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hBkg,   "Background", "l");
    leg->AddEntry(hGamma, "Gamma", "l");
    leg->Draw();

    std::string outName = OUTPUT_PATH + "gamma_spectrum_" + recLabel(gammaRecLen) + ".png";
    c1->SaveAs(outName.c_str());
    std::cout << "Saved: " << outName << std::endl;

    return 0;
}