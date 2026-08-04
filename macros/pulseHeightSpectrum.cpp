// pulseHeightSpectrum.cpp
// Compares the pulse-height spectrum of two runs (e.g. background vs DTG-on).
// Unit (ADC counts or mV) is selectable via command-line argument; defaults
// to ADC if not specified. Produces:
//   1. spectrum (rate, events/s)
//   2. counts (raw events, Background scaled to DTG-on's exposure time)
//   3. ratio (DTG-on / Background)
//   4. excess rate (DTG-on - Background, Hz)
//   5. excess counts (DTG-on - Background, raw events, Background scaled)
//
// Usage:
//   ./pulseHeightSpectrum <background_file.root> <dtgon_file.root> [unit: adc|mv]

#include "Waveform.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

static const std::string DATA_PATH   = "../data/testruns/";
static const std::string OUTPUT_PATH = "results/";

static const bool USE_LOG_Y = true;
static const int  SPECTRUM_N_BINS = 300;
static const int  RATIO_N_BINS    = 50;

static const double AMP_MIN_ADC    = 0.0;
static const double AMP_MAX_ADC    = 5000.0;
static const double AMP_CUTOFF_ADC = 700.0;

struct EventData {
    Long64_t nEntries;
    double duration;
    int recordLength;
    std::vector<double> ampsADC;
};

EventData readAllAmps(const std::string& filename)
{
    std::string fullPath = DATA_PATH + filename;
    TFile* f = TFile::Open(fullPath.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: could not open " << fullPath << std::endl;
        std::exit(1);
    }

    TTree* tree = (TTree*)f->Get("waveforms");
    if (!tree) {
        std::cerr << "ERROR: could not find tree 'waveforms' in " << fullPath << std::endl;
        std::exit(1);
    }

    Int_t recordLength;
    int waveform[20000];
    Double_t unixTime;

    tree->SetBranchAddress("recordLength", &recordLength);
    tree->SetBranchAddress("waveform", waveform);
    tree->SetBranchAddress("unixTime", &unixTime);

    Long64_t nEntries = tree->GetEntries();
    if (nEntries == 0) {
        std::cerr << "ERROR: no entries in " << fullPath << std::endl;
        std::exit(1);
    }

    tree->GetEntry(0);
    double t_start = unixTime;
    int recLen = recordLength;
    tree->GetEntry(nEntries - 1);
    double t_end = unixTime;
    double duration = t_end - t_start;

    EventData data;
    data.nEntries = nEntries;
    data.duration = duration;
    data.recordLength = recLen;
    data.ampsADC.resize(nEntries);

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, false);
        data.ampsADC[i] = w.getMaxAmpADC();
    }

    f->Close();
    std::cout << filename << ": " << nEntries << " events, duration = "
              << duration << " s, recordLength = " << recLen << std::endl;
    return data;
}

std::string recLabel(int recordLength)
{
    int recNs = std::round(recordLength * 2.0);
    int recUs = std::round(recNs / 1000.0);
    return std::to_string(recUs) + "us";
}

void makeSpectrumPlot(const EventData& bkgData, const EventData& dtgonData,
                       double ampMin, double ampMax, const std::string& unitLabel,
                       double unitScale)
{
    TH1D* hBkg = new TH1D("hBkg", "Background", SPECTRUM_N_BINS, ampMin, ampMax);
    hBkg->SetLineColor(kBlack);
    hBkg->SetLineWidth(2);
    for (double a : bkgData.ampsADC) hBkg->Fill(a * unitScale);
    hBkg->Scale(1.0 / bkgData.duration);

    TH1D* hDtgOn = new TH1D("hDtgOn", "DTG-on", SPECTRUM_N_BINS, ampMin, ampMax);
    hDtgOn->SetLineColor(kRed);
    hDtgOn->SetLineWidth(2);
    for (double a : dtgonData.ampsADC) hDtgOn->Fill(a * unitScale);
    hDtgOn->Scale(1.0 / dtgonData.duration);

    TCanvas* c1 = new TCanvas("c1", "Pulse Height Spectrum", 900, 700);
    c1->SetGrid();
    if (USE_LOG_Y) c1->SetLogy();

    double maxVal = std::max(hBkg->GetMaximum(), hDtgOn->GetMaximum());
    hBkg->SetMaximum(maxVal * 1.5);
    if (USE_LOG_Y) hBkg->SetMinimum(0.1 / std::max(bkgData.duration, dtgonData.duration));

    hBkg->GetXaxis()->SetTitle(("Max Amplitude (" + unitLabel + ")").c_str());
    hBkg->GetXaxis()->CenterTitle();
    hBkg->GetYaxis()->SetTitle("Events / second");
    hBkg->GetYaxis()->CenterTitle();

    hBkg->Draw("HIST");
    hDtgOn->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hBkg,   "Background", "l");
    leg->AddEntry(hDtgOn, "DTG-on",     "l");
    leg->Draw();

    std::string outName = OUTPUT_PATH + "spectrum_" + recLabel(bkgData.recordLength) + "_dtgandbkg.png";
    c1->SaveAs(outName.c_str());
    std::cout << "Saved: " << outName << std::endl;
}

void makeCountsPlot(const EventData& bkgData, const EventData& dtgonData,
                     double ampMin, double ampMax, const std::string& unitLabel,
                     double unitScale, double durationScale)
{
    TH1D* hBkg = new TH1D("hBkgCounts", "Background", SPECTRUM_N_BINS, ampMin, ampMax);
    hBkg->SetLineColor(kBlack);
    hBkg->SetLineWidth(2);
    for (double a : bkgData.ampsADC) hBkg->Fill(a * unitScale);
    hBkg->Scale(durationScale);

    TH1D* hDtgOn = new TH1D("hDtgOnCounts", "DTG-on", SPECTRUM_N_BINS, ampMin, ampMax);
    hDtgOn->SetLineColor(kRed);
    hDtgOn->SetLineWidth(2);
    for (double a : dtgonData.ampsADC) hDtgOn->Fill(a * unitScale);

    TCanvas* c4 = new TCanvas("c4", "Counts", 900, 700);
    c4->SetGrid();
    if (USE_LOG_Y) c4->SetLogy();

    double maxVal = std::max(hBkg->GetMaximum(), hDtgOn->GetMaximum());
    hBkg->SetMaximum(maxVal * 1.5);
    if (USE_LOG_Y) hBkg->SetMinimum(0.5);

    hBkg->GetXaxis()->SetTitle(("Max Amplitude (" + unitLabel + ")").c_str());
    hBkg->GetXaxis()->CenterTitle();
    hBkg->GetYaxis()->SetTitle("Events");
    hBkg->GetYaxis()->CenterTitle();

    hBkg->Draw("HIST");
    hDtgOn->Draw("HIST SAME");

    TLegend* leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(hBkg,   "Background", "l");
    leg->AddEntry(hDtgOn, "DTG-on",               "l");
    leg->Draw();

    std::string outName = OUTPUT_PATH + "counts_" + recLabel(bkgData.recordLength) + "_dtgandbkg.png";
    c4->SaveAs(outName.c_str());
    std::cout << "Saved: " << outName << std::endl;
}

void makeExcessCountsPlot(const EventData& bkgData, const EventData& dtgonData,
                           double ampMin, double ampMax, const std::string& unitLabel,
                           double unitScale, double durationScale)
{
    TH1D* hBkgRaw   = new TH1D("hBkgRawExC",   "Background raw", RATIO_N_BINS, ampMin, ampMax);
    TH1D* hDtgOnRaw = new TH1D("hDtgOnRawExC", "DTG-on raw",     RATIO_N_BINS, ampMin, ampMax);

    for (double a : bkgData.ampsADC)   hBkgRaw->Fill(a * unitScale);
    for (double a : dtgonData.ampsADC) hDtgOnRaw->Fill(a * unitScale);

    TH1D* hExcessCounts = new TH1D("hExcessCounts", "Excess Counts", RATIO_N_BINS, ampMin, ampMax);

    for (int bin = 1; bin <= RATIO_N_BINS; bin++) {
        double nBkgRaw   = hBkgRaw->GetBinContent(bin);
        double nDtgOnRaw = hDtgOnRaw->GetBinContent(bin);

        // Background scaled to DTG-on's exposure time, same as the counts plot.
        double nBkgScaled = nBkgRaw * durationScale;

        double excess = nDtgOnRaw - nBkgScaled;
        hExcessCounts->SetBinContent(bin, excess);

        // Error propagation: sigma(nBkgRaw)=sqrt(nBkgRaw), scaled by the same
        // durationScale factor; sigma(nDtgOnRaw)=sqrt(nDtgOnRaw). Independent,
        // so variances add.
        double sigmaBkgScaled = std::sqrt(nBkgRaw) * durationScale;
        double sigmaDtgOn     = std::sqrt(nDtgOnRaw);
        double error = std::sqrt(sigmaDtgOn*sigmaDtgOn + sigmaBkgScaled*sigmaBkgScaled);
        hExcessCounts->SetBinError(bin, error);
    }

    TCanvas* c5 = new TCanvas("c5", "Excess Counts", 900, 700);
    c5->SetGrid();

    hExcessCounts->SetLineColor(kMagenta+2);
    hExcessCounts->SetMarkerColor(kMagenta+2);
    hExcessCounts->SetMarkerStyle(20);
    hExcessCounts->SetLineWidth(2);

    hExcessCounts->GetXaxis()->SetTitle(("Max Amplitude (" + unitLabel + ")").c_str());
    hExcessCounts->GetXaxis()->CenterTitle();
    hExcessCounts->GetYaxis()->SetTitle("Excess Events: DTG-on - Background");
    hExcessCounts->GetYaxis()->CenterTitle();

    hExcessCounts->Draw("E1");

    TLine* zeroLine = new TLine(ampMin, 0.0, ampMax, 0.0);
    zeroLine->SetLineStyle(2);
    zeroLine->SetLineColor(kGray+2);
    zeroLine->Draw("SAME");

    std::string outName = OUTPUT_PATH + "excess_counts_" + recLabel(bkgData.recordLength) + "_dtgandbkg.png";
    c5->SaveAs(outName.c_str());
    std::cout << "Saved: " << outName << std::endl;
}

void makeRatioPlot(const EventData& bkgData, const EventData& dtgonData,
                    double ampMin, double ampMax, const std::string& unitLabel,
                    double unitScale)
{
    TH1D* hBkgRaw   = new TH1D("hBkgRaw",   "Background raw", RATIO_N_BINS, ampMin, ampMax);
    TH1D* hDtgOnRaw = new TH1D("hDtgOnRaw", "DTG-on raw",     RATIO_N_BINS, ampMin, ampMax);

    for (double a : bkgData.ampsADC)   hBkgRaw->Fill(a * unitScale);
    for (double a : dtgonData.ampsADC) hDtgOnRaw->Fill(a * unitScale);

    TH1D* hRatio = new TH1D("hRatio", "Ratio", RATIO_N_BINS, ampMin, ampMax);

    for (int bin = 1; bin <= RATIO_N_BINS; bin++) {
        double nBkg   = hBkgRaw->GetBinContent(bin);
        double nDtgOn = hDtgOnRaw->GetBinContent(bin);

        if (nBkg <= 0 || nDtgOn <= 0) {
            hRatio->SetBinContent(bin, 0);
            hRatio->SetBinError(bin, 0);
            continue;
        }

        double rateBkg   = nBkg   / bkgData.duration;
        double rateDtgOn = nDtgOn / dtgonData.duration;
        double ratio = rateDtgOn / rateBkg;

        double relError = std::sqrt(1.0/nDtgOn + 1.0/nBkg);
        double error = ratio * relError;

        hRatio->SetBinContent(bin, ratio);
        hRatio->SetBinError(bin, error);
    }

    TCanvas* c2 = new TCanvas("c2", "Ratio Spectrum", 900, 700);
    c2->SetGrid();

    hRatio->SetLineColor(kBlue);
    hRatio->SetMarkerColor(kBlue);
    hRatio->SetMarkerStyle(20);
    hRatio->SetLineWidth(2);

    hRatio->GetXaxis()->SetTitle(("Max Amplitude (" + unitLabel + ")").c_str());
    hRatio->GetXaxis()->CenterTitle();
    hRatio->GetYaxis()->SetTitle("DTG-on / Background");
    hRatio->GetYaxis()->CenterTitle();

    hRatio->Draw("E1");

    TLine* line = new TLine(ampMin, 1.0, ampMax, 1.0);
    line->SetLineStyle(2);
    line->SetLineColor(kGray+2);
    line->Draw("SAME");

    std::string outName = OUTPUT_PATH + "ratio_spectrum_" + recLabel(bkgData.recordLength) + "_dtgandbkg.png";
    c2->SaveAs(outName.c_str());
    std::cout << "Saved: " << outName << std::endl;
}

void makeExcessRatePlot(const EventData& bkgData, const EventData& dtgonData,
                         double ampMin, double ampMax, double ampCutoff,
                         const std::string& unitLabel, double unitScale)
{
    TH1D* hBkgRaw   = new TH1D("hBkgRawEx",   "Background raw", RATIO_N_BINS, ampMin, ampMax);
    TH1D* hDtgOnRaw = new TH1D("hDtgOnRawEx", "DTG-on raw",     RATIO_N_BINS, ampMin, ampMax);

    for (double a : bkgData.ampsADC)   hBkgRaw->Fill(a * unitScale);
    for (double a : dtgonData.ampsADC) hDtgOnRaw->Fill(a * unitScale);

    TH1D* hExcess = new TH1D("hExcess", "Excess Rate", RATIO_N_BINS, ampMin, ampMax);

    for (int bin = 1; bin <= RATIO_N_BINS; bin++) {
        double nBkg   = hBkgRaw->GetBinContent(bin);
        double nDtgOn = hDtgOnRaw->GetBinContent(bin);

        double rateBkg   = nBkg   / bkgData.duration;
        double rateDtgOn = nDtgOn / dtgonData.duration;
        double excess = rateDtgOn - rateBkg;

        double sigmaDtgOn = std::sqrt(nDtgOn) / dtgonData.duration;
        double sigmaBkg   = std::sqrt(nBkg)   / bkgData.duration;
        double error = std::sqrt(sigmaDtgOn*sigmaDtgOn + sigmaBkg*sigmaBkg);

        hExcess->SetBinContent(bin, excess);
        hExcess->SetBinError(bin, error);
    }

    TCanvas* c3 = new TCanvas("c3", "Excess Rate", 900, 700);
    c3->SetGrid();

    hExcess->SetLineColor(kGreen+2);
    hExcess->SetMarkerColor(kGreen+2);
    hExcess->SetMarkerStyle(20);
    hExcess->SetLineWidth(2);

    hExcess->GetXaxis()->SetTitle(("Max Amplitude (" + unitLabel + ")").c_str());
    hExcess->GetXaxis()->CenterTitle();
    hExcess->GetYaxis()->SetTitle("Excess Rate: DTG-on - Background (Hz)");
    hExcess->GetYaxis()->CenterTitle();

    hExcess->Draw("E1");

    TLine* zeroLine = new TLine(ampMin, 0.0, ampMax, 0.0);
    zeroLine->SetLineStyle(2);
    zeroLine->SetLineColor(kGray+2);
    zeroLine->Draw("SAME");

    TLine* cutoffLine = new TLine(ampCutoff, hExcess->GetMinimum(),
                                   ampCutoff, hExcess->GetMaximum());
    cutoffLine->SetLineStyle(2);
    cutoffLine->SetLineColor(kRed);
    cutoffLine->Draw("SAME");

    std::string outName = OUTPUT_PATH + "excess_rate_" + recLabel(bkgData.recordLength) + "_dtgandbkg.png";
    c3->SaveAs(outName.c_str());
    std::cout << "Saved: " << outName << std::endl;
}

void computeExcessRate(const EventData& bkgData, const EventData& dtgonData,
                        double ampMin, double ampMax, double ampCutoff,
                        const std::string& unitLabel, double unitScale)
{
    TH1D* hBkgRaw   = new TH1D("hBkgRawTotal",   "", RATIO_N_BINS, ampMin, ampMax);
    TH1D* hDtgOnRaw = new TH1D("hDtgOnRawTotal", "", RATIO_N_BINS, ampMin, ampMax);

    for (double a : bkgData.ampsADC)   hBkgRaw->Fill(a * unitScale);
    for (double a : dtgonData.ampsADC) hDtgOnRaw->Fill(a * unitScale);

    std::cout << "\n--- Raw bin contents, first 10 bins (" << unitLabel << ") ---" << std::endl;
    for (int bin = 1; bin <= 10; bin++) {
        double lowEdge  = hBkgRaw->GetBinLowEdge(bin);
        double highEdge = lowEdge + hBkgRaw->GetBinWidth(bin);
        double nBkg   = hBkgRaw->GetBinContent(bin);
        double nDtgOn = hDtgOnRaw->GetBinContent(bin);
        std::cout << "Bin " << bin << " [" << lowEdge << "-" << highEdge << " " << unitLabel << "]: "
                  << "nBkg=" << nBkg << "  nDtgOn=" << nDtgOn << std::endl;
    }

    double totalExcess = 0;
    double totalVariance = 0;

    for (int bin = 1; bin <= RATIO_N_BINS; bin++) {
        double binCenter = hBkgRaw->GetBinCenter(bin);
        if (binCenter <= ampCutoff) continue;

        double nBkg   = hBkgRaw->GetBinContent(bin);
        double nDtgOn = hDtgOnRaw->GetBinContent(bin);

        double rateBkg   = nBkg   / bkgData.duration;
        double rateDtgOn = nDtgOn / dtgonData.duration;

        double excessBin = rateDtgOn - rateBkg;
        totalExcess += excessBin;

        double sigmaDtgOn = std::sqrt(nDtgOn) / dtgonData.duration;
        double sigmaBkg   = std::sqrt(nBkg)   / bkgData.duration;
        totalVariance += sigmaDtgOn*sigmaDtgOn + sigmaBkg*sigmaBkg;
    }

    double totalError = std::sqrt(totalVariance);

    std::cout << "\nIntegrated excess rate above " << ampCutoff << " " << unitLabel << ": "
              << totalExcess << " +/- " << totalError << " Hz" << std::endl;

    delete hBkgRaw;
    delete hDtgOnRaw;
}

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc < 3 || argc > 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <background_file.root> <dtgon_file.root> [unit: adc|mv]" << std::endl;
        return 1;
    }

    std::string bkgFile   = argv[1];
    std::string dtgonFile = argv[2];
    std::string unit = (argc == 4) ? argv[3] : "adc";

    std::transform(unit.begin(), unit.end(), unit.begin(), ::tolower);
    if (unit != "adc" && unit != "mv") {
        std::cerr << "ERROR: unit must be 'adc' or 'mv'" << std::endl;
        return 1;
    }

    double unitScale = (unit == "mv") ? Waveform::MV_PER_COUNT : 1.0;
    std::string unitLabel = (unit == "mv") ? "mV" : "ADC";

    double ampMin    = AMP_MIN_ADC    * unitScale;
    double ampMax    = AMP_MAX_ADC    * unitScale;
    double ampCutoff = AMP_CUTOFF_ADC * unitScale;

    gSystem->mkdir(OUTPUT_PATH.c_str(), true);

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    EventData bkgData   = readAllAmps(bkgFile);
    EventData dtgonData = readAllAmps(dtgonFile);

    double durationScale = dtgonData.duration / bkgData.duration;

    makeSpectrumPlot(bkgData, dtgonData, ampMin, ampMax, unitLabel, unitScale);
    makeCountsPlot(bkgData, dtgonData, ampMin, ampMax, unitLabel, unitScale, durationScale);
    makeExcessCountsPlot(bkgData, dtgonData, ampMin, ampMax, unitLabel, unitScale, durationScale);
    makeRatioPlot(bkgData, dtgonData, ampMin, ampMax, unitLabel, unitScale);
    makeExcessRatePlot(bkgData, dtgonData, ampMin, ampMax, ampCutoff, unitLabel, unitScale);
    computeExcessRate(bkgData, dtgonData, ampMin, ampMax, ampCutoff, unitLabel, unitScale);

    return 0;
}