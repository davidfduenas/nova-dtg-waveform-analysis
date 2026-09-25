// NeutronCurveOverlay.cpp
#include "NeutronCurveOverlay.h"
#include "Waveform.h"

#include <cmath>
#include <vector>

#include <TH1D.h>
#include <TCanvas.h>
#include <TLatex.h>

NeutronCurveOverlay::NeutronCurveOverlay() {}

double NeutronCurveOverlay::evaluate(double tNs) const
{
    if (tNs < 0) return 0.0;

    double raw = W_FAST * std::exp(-tNs / T_FAST)
               + W_MID  * std::exp(-tNs / T_MID)
               + W_SLOW * std::exp(-tNs / T_SLOW);

    return raw / W_SUM; // normalized so evaluate(0) == 1.0
}

void NeutronCurveOverlay::draw(const Waveform& wf, int explicitStartSample, double nsPerSample,
                                long eventIndex, TCanvas* c, const std::string& pdfName,
                                const std::string& extraLabel) const
{
    const std::vector<double>& subtracted = wf.getSubtractedADC();
    int n = (int)subtracted.size();

    // --- Draw the waveform itself, same convention as Waveform::drawWithThreshold() ---
    std::string hname = "hNeutronCurve" + std::to_string(eventIndex);
    TH1D* h = new TH1D(hname.c_str(),
                        ("Event " + std::to_string(eventIndex) + ";Time (ns);Baseline - Signal (ADC)").c_str(),
                        n, 0, n * nsPerSample);
    for (int s = 0; s < n; s++) h->SetBinContent(s + 1, subtracted[s]);
    h->SetLineColor(kBlue);
    h->SetLineWidth(1);
    h->Draw("HIST");

    // --- Build the reference curve, anchored at explicitStartSample and
    //     scaled so its value at t=0 (i.e. at explicitStartSample) equals
    //     the event's own peak amplitude. Pure visual overlay -- no fit. ---
    std::string cname = "hRefCurve" + std::to_string(eventIndex);
    TH1D* curve = new TH1D(cname.c_str(), "", n, 0, n * nsPerSample);

    double peakAmp = wf.getMaxAmpADC();
    if (explicitStartSample >= 0 && explicitStartSample < n) {
        for (int s = 0; s < n; s++) {
            double tNs = (s - explicitStartSample) * nsPerSample;
            double val = (s < explicitStartSample) ? 0.0 : peakAmp * evaluate(tNs);
            curve->SetBinContent(s + 1, val);
        }
        curve->SetLineColor(kGreen + 2);
        curve->SetLineStyle(2);
        curve->SetLineWidth(2);
        curve->Draw("HIST SAME");
    }
    // If explicitStartSample is invalid, the curve is left empty (all
    // zeros) rather than drawn at a guessed anchor -- caller should check
    // findNeutronLikeStretch()/countSeparatePulses() succeeded before
    // calling draw().

    if (!extraLabel.empty()) {
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.035);
        latex.SetTextColor(kBlack);
        latex.DrawLatex(0.15, 0.85, extraLabel.c_str());
    }

    c->Print(pdfName.c_str());

    delete curve;
    delete h;
}