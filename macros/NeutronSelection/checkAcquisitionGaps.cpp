// checkAcquisitionGaps.cpp
// Usage: ./checkAcquisitionGaps <file.root>
// Reconstructs real elapsed time from triggerTimeStamp (31-bit counter,
// 8 ns/tick, rolls over every ~17.2 s) by unwrapping the rollovers, then
// plots reconstructed time vs event index. Flags "breaks" -- gaps between
// consecutive events much longer than the typical spacing -- which would
// indicate the board was busy (buffer full / readout stall) and not
// listening for triggers, rather than a real absence of pulses.
//
// Only reads triggerTimeStamp (not the waveform array), so this runs fast
// even on multi-million-event files.

#include "Utilities.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TError.h>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

static const double NS_PER_TICK = 8.0;
static const Long64_t COUNTER_RANGE = 1LL << 31; // 31-bit counter

// A gap is flagged as a "break" if it's this many times longer than the
// typical (median) spacing between consecutive events.
static const double JUMP_MULT = 20.0;

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

    // Only read triggerTimeStamp -- skip the (huge) waveform branch entirely.
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("triggerTimeStamp", 1);

    Long64_t rawTick;
    tree->SetBranchAddress("triggerTimeStamp", &rawTick);

    Long64_t nEntries = tree->GetEntries();
    std::cout << inFile << ": " << nEntries << " total events" << std::endl;

    // --- Unwrap the 31-bit rollover ---
    std::vector<double> timeSec(nEntries);
    Long64_t prevRaw = 0;
    Long64_t offset = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        if (i > 0 && rawTick < prevRaw - (COUNTER_RANGE / 2)) {
            offset += COUNTER_RANGE; // genuine forward rollover
        }
        Long64_t unwrapped = rawTick + offset;
        timeSec[i] = unwrapped * NS_PER_TICK / 1e9;
        prevRaw = rawTick;
    }

    // --- Compute inter-event gaps and the median (typical) gap ---
    std::vector<double> gaps(nEntries - 1);
    for (Long64_t i = 1; i < nEntries; i++) gaps[i - 1] = timeSec[i] - timeSec[i - 1];

    std::vector<double> sortedGaps = gaps;
    std::sort(sortedGaps.begin(), sortedGaps.end());
    double medianGap = sortedGaps[sortedGaps.size() / 2];

    // --- Find breaks and build segments ---
    std::vector<Long64_t> segmentStarts = {0};
    for (Long64_t i = 1; i < nEntries; i++) {
        if (gaps[i - 1] > JUMP_MULT * medianGap) segmentStarts.push_back(i);
    }
    segmentStarts.push_back(nEntries); // sentinel end

    int nSegments = (int)segmentStarts.size() - 1;
    double totalLiveTime = 0;
    long totalBreakEvents = 0;

    std::vector<std::string> segmentLines;
    for (int s = 0; s < nSegments; s++) {
        Long64_t start = segmentStarts[s];
        Long64_t end = segmentStarts[s + 1] - 1;
        double dt = timeSec[end] - timeSec[start];
        long nEv = end - start + 1;
        double rate = (dt > 0) ? nEv / dt : 0;
        totalLiveTime += dt;

        char buf[128];
        snprintf(buf, sizeof(buf), "Seg %d: rate %.2f Hz | dt %.3f s | ev[%lld..%lld]",
                  s, rate, dt, start, end);
        segmentLines.push_back(buf);
    }

    double weightedAvgRate = (totalLiveTime > 0) ? nEntries / totalLiveTime : 0;

    std::cout << "\nSegments (breaks > " << JUMP_MULT << "x median gap of "
              << medianGap * 1e6 << " us): " << nSegments << std::endl;
    for (auto& line : segmentLines) std::cout << "  " << line << std::endl;
    std::cout << "Weighted avg rate: " << weightedAvgRate << " Hz (live time "
              << totalLiveTime << " s)" << std::endl;

    // --- Plot: reconstructed time vs event index ---
    std::string outDir = Utilities::makeOutputDir(argv[0], argv[1]);
    gSystem->mkdir(outDir.c_str(), true);
    std::string pdfName = outDir + "acquisitionGaps.pdf";

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(1);

    std::vector<double> eventIdx(nEntries);
    for (Long64_t i = 0; i < nEntries; i++) eventIdx[i] = (double)i;

    TGraph* g = new TGraph(nEntries, eventIdx.data(), timeSec.data());
    g->SetTitle((inFile + ";Event index;Reconstructed time (s)").c_str());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.3);
    g->SetMarkerColor(kBlue);
    g->SetLineColor(kBlue);

    TCanvas* c = new TCanvas("c", "", 1000, 700);
    c->SetGrid();
    c->Print((pdfName + "[").c_str());
    g->Draw("AP");

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.025);
    latex.SetTextColor(kBlack);
    double y = 0.85;
    char header[128];
    snprintf(header, sizeof(header), "jump_mult = %.1f | segments = %d", JUMP_MULT, nSegments);
    latex.DrawLatex(0.15, y, header); y -= 0.035;
    snprintf(header, sizeof(header), "Weighted avg rate = %.2f Hz (live time %.2f s)", weightedAvgRate, totalLiveTime);
    latex.DrawLatex(0.15, y, header); y -= 0.045;

    int maxLinesToShow = 15; // avoid running off the canvas for many-segment files
    for (int i = 0; i < (int)segmentLines.size() && i < maxLinesToShow; i++) {
        latex.DrawLatex(0.15, y, segmentLines[i].c_str());
        y -= 0.03;
    }
    if ((int)segmentLines.size() > maxLinesToShow) {
        latex.DrawLatex(0.15, y, Form("... and %d more segments (see console output)",
                                      (int)segmentLines.size() - maxLinesToShow));
    }

    c->Print(pdfName.c_str());

    // --- Second plot: gap size itself vs event index (not cumulative time).
    // A normal event sits near zero regardless of total run length; any
    // real stall shows up as an obvious spike, unlike the cumulative-time
    // plot where a long run can hide small-looking (but real) gaps. ---
    double minGap = *std::min_element(gaps.begin(), gaps.end());
    double maxGap = *std::max_element(gaps.begin(), gaps.end());
    if (minGap <= 0) minGap = 1e-7; // guard against log(0)

    int nLogBins = 60;
    std::vector<double> binEdges(nLogBins + 1);
    double logMin = std::log10(minGap);
    double logMax = std::log10(maxGap);
    for (int i = 0; i <= nLogBins; i++) {
        binEdges[i] = std::pow(10.0, logMin + (logMax - logMin) * i / nLogBins);
    }

    TH1D* hGap = new TH1D("hGap", (inFile + " -- gap duration distribution;Gap since previous event (s);Count").c_str(),
                           nLogBins, binEdges.data());
    for (double g : gaps) hGap->Fill(std::max(g, minGap));

    c->Clear();
    c->SetLogy();
    c->SetLogx();
    hGap->SetLineColor(kRed + 1);
    hGap->SetLineWidth(2);
    hGap->Draw("HIST");

    TLatex latex2;
    latex2.SetNDC();
    latex2.SetTextSize(0.03);
    latex2.SetTextColor(kBlack);
    latex2.DrawLatex(0.15, 0.85, Form("Median gap = %.3g s | break threshold = %.3g s (%.1fx median)",
                                       medianGap, JUMP_MULT * medianGap, JUMP_MULT));

    c->Print(pdfName.c_str());
    c->Print((pdfName + "]").c_str());

    delete hGap;
    std::cout << "Saved: " << pdfName << std::endl;

    delete g;
    f->Close();
    return 0;
}