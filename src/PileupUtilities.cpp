// PileupUtilities.cpp
#include "PileupUtilities.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>

#include <algorithm>
#include <cstdio>
#include <iostream>

PileupUtilities::PileupUtilities(const std::string& outDir, const Settings& settings)
    : dir(outDir), cfg(settings), nDropped(0), guardNs(0.0),
      countByHeight(settings.tableEdges.size() + 1, 0)
{
    hHeight = new TH1D("hPileupJumpHeight", "Jumps after the first pulse;Jump height (ADC);Jumps",
                       cfg.heightBins, 0, cfg.heightMax);
    hRatio  = new TH1D("hPileupJumpRatio", "Jump height relative to the first pulse;Height / first pulse height;Jumps",
                       cfg.ratioBins, 0, cfg.ratioMax);

    char title[200];
    hFirstLater = new TH2D("hPileupFirstLater",
                           "First climb against the biggest later climb in the same event;First climb height (ADC);Biggest later climb height (ADC)",
                           cfg.firstBins, 0, cfg.firstMax, cfg.laterBins, 0, cfg.laterMax);
    hLaterMax = new TH1D("hPileupLaterMax", "", cfg.laterBins, 0, cfg.laterMax);
    snprintf(title, sizeof(title), "Time from the first pulse to later pulses;Gap (ns);Pulses");
    hGap = new TH1D("hPileupGap", title, cfg.gapBins, 0, cfg.gapMax);
    hFirstAll = new TH1D("hPileupFirstAll", "", cfg.firstBins, 0, cfg.firstMax);
    hFirstWithLater = new TH1D("hPileupFirstWithLater", "", cfg.firstBins, 0, cfg.firstMax);
    hFracNoise = new TH2D("hPileupFracNoise",
                          "Every climb after the first pulse;Height as a fraction of the event's peak;Height divided by the noise before it",
                          cfg.fracBins, 0, cfg.fracMax, cfg.sigmaBins, 0, cfg.sigmaMax);

    // Keep the histograms out of the input file, which is closed before this object is destroyed.
    hHeight->SetDirectory(nullptr);
    hRatio->SetDirectory(nullptr);
    hFirstLater->SetDirectory(nullptr);
    hLaterMax->SetDirectory(nullptr);
    hGap->SetDirectory(nullptr);
    hFirstAll->SetDirectory(nullptr);
    hFirstWithLater->SetDirectory(nullptr);
    hFracNoise->SetDirectory(nullptr);
}

PileupUtilities::~PileupUtilities()
{
    delete hHeight;
    delete hRatio;
    delete hFirstLater;
    delete hLaterMax;
    delete hGap;
    delete hFirstAll;
    delete hFirstWithLater;
    delete hFracNoise;
}

void PileupUtilities::add(long eventIndex, const std::vector<double>& shape,
                          const std::vector<PileupFinder::Pulse>& pulses, double psd)
{
    if ((long)events.size() >= cfg.maxStored) { nDropped++; return; }
    Stored s;
    s.index = eventIndex;
    s.psd = psd;
    s.pulses = pulses;
    s.shape.assign(shape.begin(), shape.end());
    events.push_back(s);
}

void PileupUtilities::write()
{
    TCanvas* c = new TCanvas("cPileup", "", cfg.canvasWidth, cfg.canvasHeight);
    c->SetGrid();
    writeCategory(c, "1pulse", 1, 1);
    writeCategory(c, "2pulses", 2, 2);
    writeCategory(c, "3pluspulses", 3, 1000000);
    delete c;
}

void PileupUtilities::writeCategory(TCanvas* c, const std::string& tag, int minP, int maxP)
{
    int part = 0, pagesInPart = 0;
    long total = 0;
    std::string path;

    for (size_t k = 0; k < events.size(); k++) {
        int np = (int)events[k].pulses.size();
        if (np < minP || np > maxP) continue;

        if (pagesInPart == 0) {
            part++;
            char num[16];
            snprintf(num, sizeof(num), "%03d", part);
            path = dir + "pileup_" + tag + "_part" + num + ".pdf";
            c->Print((path + "[").c_str());
        }

        drawEvent(c, events[k], path);
        pagesInPart++;
        total++;

        if (pagesInPart >= cfg.pagesPerPdf) {
            c->Print((path + "]").c_str());
            pagesInPart = 0;
        }
    }
    if (pagesInPart > 0) c->Print((path + "]").c_str());

    std::cout << tag << ": " << total << " events, " << part << " PDF file(s) in " << dir << std::endl;
}

void PileupUtilities::drawEvent(TCanvas* c, const Stored& s, const std::string& pdfName)
{
    int n = (int)s.shape.size();

    char title[128];
    snprintf(title, sizeof(title), "Event %ld, %d %s;Time (ns);Baseline - Signal (ADC)",
             s.index, (int)s.pulses.size(), s.pulses.size() == 1 ? "pulse" : "pulses");

    TH1D* h = new TH1D(("hPileup" + std::to_string(s.index)).c_str(), title, n, 0, n * cfg.nsPerSample);
    for (int j = 0; j < n; j++) h->SetBinContent(j + 1, s.shape[j]);
    h->SetStats(0);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h->SetMaximum(h->GetMaximum() * cfg.headroom); // headroom for the text
    h->Draw("HIST");

    double yMin = h->GetMinimum();
    double yMax = h->GetMaximum();

    std::vector<TLine*> lines;
    for (size_t k = 0; k < s.pulses.size(); k++) {
        double x = s.pulses[k].onsetSample * cfg.nsPerSample;
        TLine* l = new TLine(x, yMin, x, yMax);
        l->SetLineColor(kRed);
        l->SetLineStyle(2);
        l->SetLineWidth(1);
        l->Draw();
        lines.push_back(l);
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.03);
    latex.SetTextColor(kBlack);

    char buf[160];
    snprintf(buf, sizeof(buf), "PSD = %.3f", s.psd);
    latex.DrawLatex(0.15, 0.86, buf);

    double y = 0.86;
    for (size_t k = 0; k < s.pulses.size() && (int)k < cfg.maxPulsesListed; k++) {
        snprintf(buf, sizeof(buf), "pulse %d: t = %.0f ns, height = %.0f ADC",
                 (int)k + 1, s.pulses[k].onsetSample * cfg.nsPerSample, s.pulses[k].height);
        latex.DrawLatex(0.50, y, buf);
        y -= 0.04;
    }

    c->Print(pdfName.c_str());

    for (size_t k = 0; k < lines.size(); k++) delete lines[k];
    delete h;
}

int PileupUtilities::heightBin(double height) const
{
    for (size_t i = 0; i < cfg.tableEdges.size(); i++) {
        if (height < cfg.tableEdges[i]) return (int)i;
    }
    return (int)cfg.tableEdges.size();
}

void PileupUtilities::addEvent(const std::vector<PileupFinder::Pulse>& wideClimbs,
                               const std::vector<PileupFinder::Pulse>& pulses,
                               double eventPeak, int guardSamples)
{
    if (wideClimbs.empty() || pulses.empty()) return;

    const PileupFinder::Pulse& first = wideClimbs[0];
    if (first.height <= 0.0) return;

    guardNs = guardSamples * cfg.nsPerSample;

    // Every climb after the first one, whether or not it counts as a pulse.
    double biggestLater = 0.0;
    for (size_t k = 1; k < wideClimbs.size(); k++) {
        const PileupFinder::Pulse& p = wideClimbs[k];
        if (p.onsetSample < first.onsetSample + guardSamples) continue; // jagged top of the first pulse

        hHeight->Fill(p.height);
        hRatio->Fill(p.height / first.height);
        countByHeight[heightBin(p.height)]++;
        if (eventPeak > 0.0 && p.noise > 0.0) hFracNoise->Fill(p.height / eventPeak, p.height / p.noise);

        if (p.height > biggestLater) biggestLater = p.height;
    }
    if (biggestLater > 0.0) {
        hLaterMax->Fill(biggestLater);
        hFirstLater->Fill(first.height, biggestLater);
    }

    // The pulses the finder counted.
    hFirstAll->Fill(pulses[0].height);
    if (pulses.size() >= 2) {
        hFirstWithLater->Fill(pulses[0].height);
        for (size_t k = 1; k < pulses.size(); k++) {
            hGap->Fill((pulses[k].onsetSample - pulses[0].onsetSample) * cfg.nsPerSample);
        }
    }
}

void PileupUtilities::writeHeightPlots(const std::string& pdfPath)
{
    TCanvas* c = new TCanvas("cPileupHeights", "", cfg.canvasWidth, cfg.canvasHeight);
    c->SetGrid();
    c->SetLogy(1);
    c->Print((pdfPath + "[").c_str());

    hHeight->SetStats(0);
    hHeight->SetLineColor(kBlack);
    hHeight->SetLineWidth(2);
    hHeight->SetMinimum(0.5);
    hHeight->Draw("HIST");
    c->Print(pdfPath.c_str());

    hRatio->SetStats(0);
    hRatio->SetLineColor(kBlack);
    hRatio->SetLineWidth(2);
    hRatio->SetMinimum(0.5);
    hRatio->Draw("HIST");
    c->Print(pdfPath.c_str());

    c->Print((pdfPath + "]").c_str());
    c->Clear();
    delete c;
}

void PileupUtilities::printHeightSummary() const
{
    std::cout << "Jumps after the first pulse (wide-open search):" << std::endl;
    std::cout << "  height (ADC)       jumps" << std::endl;
    for (size_t i = 0; i < countByHeight.size(); i++) {
        char label[48];
        if (cfg.tableEdges.empty())            snprintf(label, sizeof(label), "all");
        else if (i == 0)                       snprintf(label, sizeof(label), "below %.0f", cfg.tableEdges[0]);
        else if (i == cfg.tableEdges.size())   snprintf(label, sizeof(label), "%.0f and up", cfg.tableEdges[i - 1]);
        else                                   snprintf(label, sizeof(label), "%.0f - %.0f", cfg.tableEdges[i - 1], cfg.tableEdges[i]);
        std::cout << "  " << label << "    " << countByHeight[i] << std::endl;
    }
}

void PileupUtilities::writeStatisticsPlots(const std::string& pdfPath)
{
    TCanvas* c = new TCanvas("cPileupStats", "", cfg.canvasWidth, cfg.canvasHeight);
    c->SetGrid();
    c->Print((pdfPath + "[").c_str());

    // 1. First climb against the biggest later climb.
    c->SetLogy(0);
    c->SetLogz(1);
    hFirstLater->SetStats(0);
    hFirstLater->Draw("COLZ");
    c->Print(pdfPath.c_str());

    // 2. Events with a later climb of at least a given height (cumulative from the top).
    TH1D* hCum = new TH1D("hPileupCumulative",
                          "Events with a later climb of at least this height;Later climb height (ADC);Events",
                          cfg.laterBins, 0, cfg.laterMax);
    hCum->SetDirectory(nullptr);
    double running = hLaterMax->GetBinContent(cfg.laterBins + 1); // overflow
    for (int b = cfg.laterBins; b >= 1; b--) {
        running += hLaterMax->GetBinContent(b);
        hCum->SetBinContent(b, running);
    }
    c->SetLogz(0);
    c->SetLogy(1);
    hCum->SetStats(0);
    hCum->SetLineColor(kBlack);
    hCum->SetLineWidth(2);
    hCum->SetMinimum(0.5);
    hCum->Draw("HIST");
    c->Print(pdfPath.c_str());

    // 3. Gap between the first climb and later big climbs, with a flat reference.
    c->SetLogy(0);
    hGap->SetStats(0);
    hGap->SetLineColor(kBlack);
    hGap->SetLineWidth(2);
    hGap->SetMinimum(0);
    hGap->Draw("HIST");

    double binWidth = cfg.gapMax / cfg.gapBins;
    double binsInRange = (cfg.gapFlatEnd - guardNs) / binWidth;
    TLine* flat = nullptr;
    if (binsInRange > 0) {
        double level = hGap->GetEntries() / binsInRange;
        flat = new TLine(guardNs, level, cfg.gapFlatEnd, level);
        flat->SetLineColor(kRed);
        flat->SetLineStyle(2);
        flat->SetLineWidth(2);
        flat->Draw();
    }
    c->Print(pdfPath.c_str());

    // 4. Fraction of events with a later pulse, against the height of the first pulse.
    TH1D* hFrac = (TH1D*)hFirstWithLater->Clone("hPileupFraction");
    hFrac->SetDirectory(nullptr);
    hFrac->Divide(hFirstWithLater, hFirstAll, 1, 1, "B");
    hFrac->SetTitle("Fraction of events with a later pulse;First pulse height (ADC);Fraction of events");
    hFrac->SetStats(0);
    hFrac->SetMinimum(0);
    hFrac->SetMarkerStyle(20);
    hFrac->SetLineColor(kBlack);
    hFrac->Draw("E");

    double all = hFirstAll->GetEntries();
    TLine* mean = nullptr;
    if (all > 0) {
        double overall = hFirstWithLater->GetEntries() / all;
        mean = new TLine(0, overall, cfg.firstMax, overall);
        mean->SetLineColor(kRed);
        mean->SetLineStyle(2);
        mean->SetLineWidth(2);
        mean->Draw();
    }
    c->Print(pdfPath.c_str());

    // 5. Every climb after the first pulse: height as a fraction of the event's
    // peak against height divided by the noise before it. A later climb counts
    // as a pulse if it is right of the vertical line AND above the horizontal one.
    c->SetLogy(0);
    c->SetLogz(1);
    hFracNoise->SetStats(0);
    hFracNoise->Draw("COLZ");
    TLine* cutX = new TLine(cfg.cutFraction, 0, cfg.cutFraction, cfg.sigmaMax);
    TLine* cutY = new TLine(0, cfg.cutSigmas, cfg.fracMax, cfg.cutSigmas);
    cutX->SetLineColor(kRed);
    cutX->SetLineStyle(2);
    cutX->SetLineWidth(2);
    cutX->Draw();
    cutY->SetLineColor(kRed);
    cutY->SetLineStyle(2);
    cutY->SetLineWidth(2);
    cutY->Draw();
    c->Print(pdfPath.c_str());
    c->SetLogz(0);

    c->Print((pdfPath + "]").c_str());

    delete cutX;
    delete cutY;
    delete flat;
    delete mean;
    delete hCum;
    delete hFrac;
    c->Clear();
    delete c;
}
