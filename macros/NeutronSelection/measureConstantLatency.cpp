// measureConstantLatency.cpp
// Usage: ./measureConstantLatency <background_file.root>
//
// Example, using our actual 2us setup:
//   recordLength = 1030 (rounded up from requested 1024)
//   PostTriggerValue = 103 (read from register 0x8114)
//   N = 8
//   intendedNpost = 103 * 8 = 824
//   intendedNpre  = 1030 - 824 = 206   (this is the 20% pretrigger we asked for)
//
//   Then we measure, from real data, the sample where each event's signal
//   crosses 50 ADC (the real trigger threshold). Say the average of that
//   comes out to 132.
//
//   ConstantLatency = intendedNpre - measured = 206 - 132 = 74
//
//   So the real pretrigger is only 132 samples (~264 ns, ~12.8%), not the
//   206 samples (20%) we asked for -- the other 74 samples get eaten by
//   the board's trigger-processing delay before recording actually starts
//   counting post-trigger samples.

#include "Waveform.h"

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
#include <algorithm>
#include <cmath>

static const std::string OUTPUT_PATH = "results/";
static const std::string INPUT_DIR   = "/Users/david/DTGAnalysis/data/testruns/";

static const double POST_TRIGGER_VALUE = 103.0; // from register readback
static const double N_COEFFICIENT      = 8.0;   // 730 family
static const double THRESHOLD_ADC      = 50.0;  // hardware self-trigger threshold

int main(int argc, char** argv)
{
    gROOT->SetBatch(kTRUE);

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <background_file.root>" << std::endl;
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

    Long64_t nEntries = tree->GetEntries();
    std::cout << inFile << ": " << nEntries << " total events" << std::endl;

    std::vector<double> crossingSamples;
    crossingSamples.reserve(nEntries);
    long noCrossingCount = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        Waveform w(waveform, recordLength, true);
        const std::vector<double>& shape = w.getSubtractedADC();

        int crossingSample = -1;
        for (int s = 0; s < recordLength; s++) {
            if (shape[s] > THRESHOLD_ADC) { crossingSample = s; break; }
        }

        if (crossingSample >= 0) crossingSamples.push_back(crossingSample);
        else noCrossingCount++;
    }

    if (crossingSamples.empty()) {
        std::cerr << "No crossings found." << std::endl;
        f->Close();
        return 1;
    }

    double sum = 0;
    for (double v : crossingSamples) sum += v;
    double mean = sum / crossingSamples.size();

    double sumSq = 0;
    for (double v : crossingSamples) sumSq += (v - mean) * (v - mean);
    double stdDev = std::sqrt(sumSq / (crossingSamples.size() - 1));
    double stdErrMean = stdDev / std::sqrt((double)crossingSamples.size());

    std::sort(crossingSamples.begin(), crossingSamples.end());
    double median = crossingSamples[crossingSamples.size() / 2];

    std::cout << "\nCrossings: " << crossingSamples.size() << " / " << nEntries
              << " (no crossing: " << noCrossingCount << ")" << std::endl;
    std::cout << "mean = " << mean << " +/- " << stdErrMean << std::endl;
    std::cout << "median = " << median << std::endl;
    std::cout << "stdDev = " << stdDev << std::endl;

    double intendedNpost = POST_TRIGGER_VALUE * N_COEFFICIENT;
    double intendedNpre  = recordLength - intendedNpost;
    double constantLatency = intendedNpre - mean;

    const double NS_PER_SAMPLE = 2.0; // 500 MS/s digitizer

    std::cout << "\nintendedNpost = " << intendedNpost << std::endl;
    std::cout << "intendedNpre = " << intendedNpre << std::endl;
    std::cout << "measured Npre = " << mean << " +/- " << stdErrMean << std::endl;
    std::cout << "ConstantLatency = " << constantLatency << " +/- " << stdErrMean << std::endl;

    double measuredNpost = recordLength - mean;
    std::cout << "\nrecordLength = " << recordLength
              << " samples = " << recordLength * NS_PER_SAMPLE << " ns" << std::endl;
    std::cout << "measured pretrigger  = " << mean << " +/- " << stdErrMean
              << " samples = " << mean * NS_PER_SAMPLE << " +/- " << stdErrMean * NS_PER_SAMPLE << " ns" << std::endl;
    std::cout << "measured post-trigger = " << measuredNpost << " +/- " << stdErrMean
              << " samples = " << measuredNpost * NS_PER_SAMPLE << " +/- " << stdErrMean * NS_PER_SAMPLE << " ns" << std::endl;

    gSystem->mkdir(OUTPUT_PATH.c_str(), true);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1111);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    // Crossing sample is a plain integer, so use exactly 1 sample per bin.
    // Zoom to +/- 6 stdDev around the mean -- the real population lives in
    // a tight window; a handful of stray early/late crossings shouldn't
    // stretch the whole axis out and squash everything else flat.
    int zoomMin = (int)std::floor(mean - 6 * stdDev);
    int zoomMax = (int)std::ceil(mean + 6 * stdDev);
    if (zoomMin < 0) zoomMin = 0;
    if (zoomMax > recordLength) zoomMax = recordLength;

    long outsideZoomCount = 0;
    for (double v : crossingSamples) {
        if (v < zoomMin || v > zoomMax) outsideZoomCount++;
    }

    int nBins = zoomMax - zoomMin + 1;
    TH1D* h = new TH1D("hCrossingSample", ";Sample index;Events",
                        nBins, zoomMin - 0.5, zoomMax + 0.5);
    h->StatOverflows(kTRUE); // so the stat box matches the printed values above, not just the visible range
    for (double v : crossingSamples) h->Fill(v);

    TCanvas* c = new TCanvas("c", "", 900, 700);
    h->SetLineColor(kBlue);
    h->SetLineWidth(2);
    h->Draw("HIST");

    std::cout << "\nEvents outside zoomed plot range [" << zoomMin << ", " << zoomMax
              << "]: " << outsideZoomCount << " / " << crossingSamples.size() << std::endl;

    std::string outName = argv[1];
    for (char& ch : outName) if (ch == '.') ch = '_';
    std::string pdfName = OUTPUT_PATH + "crossingSampleDistribution_" + outName + ".pdf";
    c->Print(pdfName.c_str());
    std::cout << "Saved: " << pdfName << std::endl;

    delete h;
    f->Close();
    return 0;
}