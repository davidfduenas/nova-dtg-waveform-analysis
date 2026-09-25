// PSD.cpp
#include "PSD.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLine.h>
#include <TLatex.h>

#include <algorithm>

PSD::PSD(double postTriggerPercent_,
         double nCoefficient_,
         double constantLatencySamples_,
         double shortGateNs_,
         double longGateNs_,
         double nsPerSample_,
         double amplitudeCutoff_,
         double startShiftNs_)
    : postTriggerPercent(postTriggerPercent_),
      nCoefficient(nCoefficient_),
      constantLatencySamples(constantLatencySamples_),
      shortGateNs(shortGateNs_),
      longGateNs(longGateNs_),
      nsPerSample(nsPerSample_),
      amplitudeCutoff(amplitudeCutoff_),
      startShiftNs(startShiftNs_),
      totalCount(0),
      candidateCount(0)
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1111);
    gStyle->SetPadGridX(true);
    gStyle->SetPadGridY(true);

    hPSD = new TH1D("hPSD", ";PSD ratio;Events", 200, -0.2, 1.2);
    hPSD->SetDirectory(nullptr); // don't let the currently-open input file own this

    hPSDvsAmp = new TH2D("hPSDvsAmp", ";Amplitude (ADC);PSD ratio",
                          500, 0, 5000, 350, -0.2, 1.2);
    hPSDvsAmp->SetDirectory(nullptr); // same reason
}

PSD::~PSD()
{
    delete hPSD;
    delete hPSDvsAmp;
}

int PSD::getPulseStartSample(int recordLength) const
{
    // Npost = PostTriggerValue * N + ConstantLatency, and
    // PostTriggerValue*N == (postTriggerPercent/100) * recordLength for
    // how CAEN_DGTZ_SetPostTriggerSize converts a percentage into the
    // register (confirmed: PostTriggerValue=103 for recordLength=1030,
    // post=80 -> 0.80*1030/8 = 103 exactly).
    //
    // NOTE: the pre-gate shift (startShiftNs) is NOT applied here -- it's
    // applied once, uniformly, inside analyzeFromStart() below, so both
    // this fixed-formula anchor (analyze()) and an explicit detected
    // anchor (analyzeAt()) get the same pre-gate treatment consistently.
    double intendedNpost = (postTriggerPercent / 100.0) * recordLength;
    double npre = recordLength - intendedNpost - constantLatencySamples;
    return (int)(npre + 0.5); // round to nearest sample
}

PSD::Result PSD::analyzeFromStart(const std::vector<double>& shape, int recordLength, int startSample) const
{
    Result r;

    // Pre-gate: shift the anchor earlier by startShiftNs, applied here
    // uniformly regardless of whether startSample came from the fixed
    // formula (analyze()) or an explicit detected location (analyzeAt()).
    // Matches the reference paper's description of an explicit pre-gate
    // before the real charge integration -- captures a bit of genuine
    // baseline before the detected rise, rather than anchoring exactly at
    // the first threshold crossing (which can land almost on the peak for
    // a fast-rising gamma pulse).
    int preGateSamples = (int)(startShiftNs / nsPerSample + 0.5);
    r.pulseStartSample = startSample - preGateSamples;

    int shortGateSamples = (int)(shortGateNs / nsPerSample + 0.5);
    int longGateSamples  = (int)(longGateNs  / nsPerSample + 0.5);

    r.shortGateEndSample = r.pulseStartSample + shortGateSamples;
    r.longGateEndSample  = r.pulseStartSample + longGateSamples;

    // Clamp to the actual record so a pulse near the end of a short
    // record doesn't read past the array.
    int shortEnd = std::min(r.shortGateEndSample, recordLength);
    int longEnd  = std::min(r.longGateEndSample, recordLength);
    int start    = std::max(r.pulseStartSample, 0);

    double qShort = 0, qLong = 0;
    for (int s = start; s < longEnd; s++) {
        qLong += shape[s];
        if (s < shortEnd) qShort += shape[s];
    }

    r.qShort = qShort;
    r.qLong  = qLong;
    r.psdRatio = (qLong != 0) ? (qLong - qShort) / qLong : 0;

    return r;
}

PSD::Result PSD::analyze(const std::vector<double>& shape, int recordLength) const
{
    return analyzeFromStart(shape, recordLength, getPulseStartSample(recordLength));
}

PSD::Result PSD::analyzeAt(const std::vector<double>& shape, int recordLength, int explicitStartSample) const
{
    return analyzeFromStart(shape, recordLength, explicitStartSample);
}

void PSD::addEvent(const std::vector<double>& shape, int recordLength, double maxAmp)
{
    totalCount++;
    if (maxAmp <= amplitudeCutoff) return;
    candidateCount++;

    Result r = analyze(shape, recordLength);
    hPSD->Fill(r.psdRatio);
    hPSDvsAmp->Fill(maxAmp, r.psdRatio);
}

void PSD::addEventAt(const std::vector<double>& shape, int recordLength,
                      double maxAmp, int explicitStartSample)
{
    totalCount++;
    if (maxAmp <= amplitudeCutoff) return;
    candidateCount++;

    Result r = analyzeAt(shape, recordLength, explicitStartSample);
    hPSD->Fill(r.psdRatio);
    hPSDvsAmp->Fill(maxAmp, r.psdRatio);
}

void PSD::savePlots(const std::string& outputPath, const std::string& label) const
{
    gSystem->mkdir(outputPath.c_str(), true);

    TCanvas c1("c1", "", 900, 700);
    hPSD->Draw("HIST");
    c1.Print((outputPath + "psdRatio_" + label + ".pdf").c_str());

    TCanvas c2("c2", "", 900, 700);
    c2.SetRightMargin(0.15);
    hPSDvsAmp->Draw("COLZ");
    c2.Print((outputPath + "psdVsAmplitude_" + label + ".pdf").c_str());
}

void PSD::saveCombinedPlot(const PSD& denseSample, const PSD& sparseSample,
                            const std::string& outputPath, const std::string& label)
{
    gSystem->mkdir(outputPath.c_str(), true);

    TH2D* combined = (TH2D*)denseSample.hPSDvsAmp->Clone("hPSDvsAmpCombined");
    combined->Add(sparseSample.hPSDvsAmp);
    combined->SetStats(0);

    TCanvas c("cCombined", "", 900, 700);
    c.SetRightMargin(0.15);
    combined->Draw("COLZ");

    c.Print((outputPath + "psdVsAmplitude_combined_" + label + ".pdf").c_str());

    delete combined;
}

void PSD::drawFromStart(const std::vector<double>& shape, int recordLength, int startSample,
                         long eventIndex, TCanvas* c, const std::string& pdfName) const
{
    Result r = analyzeFromStart(shape, recordLength, startSample);

    std::string hname = "hEvt" + std::to_string(eventIndex);
    TH1D* h = new TH1D(hname.c_str(), ("Event " + std::to_string(eventIndex) + ";Time (ns);Baseline - Signal (ADC)").c_str(),
                        recordLength, 0, recordLength * nsPerSample);
    for (int s = 0; s < recordLength; s++) h->SetBinContent(s + 1, shape[s]);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h->Draw("HIST");

    double yMin = h->GetMinimum();
    double yMax = h->GetMaximum();

    double xStart = r.pulseStartSample * nsPerSample;
    double xShort = r.shortGateEndSample * nsPerSample;
    double xLong  = r.longGateEndSample * nsPerSample;

    TLine* lStart = new TLine(xStart, yMin, xStart, yMax);
    lStart->SetLineColor(kGreen + 2);
    lStart->SetLineStyle(2);
    lStart->SetLineWidth(1);
    lStart->Draw();

    TLine* lShort = new TLine(xShort, yMin, xShort, yMax);
    lShort->SetLineColor(kOrange + 1);
    lShort->SetLineStyle(2);
    lShort->SetLineWidth(1);
    lShort->Draw();

    TLine* lLong = new TLine(xLong, yMin, xLong, yMax);
    lLong->SetLineColor(kRed);
    lLong->SetLineStyle(2);
    lLong->SetLineWidth(1);
    lLong->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.035);
    latex.SetTextColor(kBlack);
    latex.DrawLatex(0.15, 0.85, Form("PSD = %.3f", r.psdRatio));
    latex.DrawLatex(0.15, 0.80, Form("Q_{short} = %.1f, Q_{long} = %.1f", r.qShort, r.qLong));

    c->Print(pdfName.c_str());

    delete lStart;
    delete lShort;
    delete lLong;
    delete h;
}

void PSD::drawEventWaveform(const std::vector<double>& shape, int recordLength,
                             long eventIndex, TCanvas* c, const std::string& pdfName) const
{
    drawFromStart(shape, recordLength, getPulseStartSample(recordLength), eventIndex, c, pdfName);
}

void PSD::drawEventWaveformAt(const std::vector<double>& shape, int recordLength,
                               int explicitStartSample, long eventIndex,
                               TCanvas* c, const std::string& pdfName) const
{
    drawFromStart(shape, recordLength, explicitStartSample, eventIndex, c, pdfName);
}

long PSD::getCandidateCount() const { return candidateCount; }
long PSD::getTotalCount()     const { return totalCount; }