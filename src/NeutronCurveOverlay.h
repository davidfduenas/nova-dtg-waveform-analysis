// NeutronCurveOverlay.h
#ifndef NEUTRONCURVEOVERLAY_H
#define NEUTRONCURVEOVERLAY_H

#include <string>

class Waveform;
class TCanvas;

// Overlays the paper's reference neutron decay curve on top of a single
// event's baseline-subtracted waveform, for visual (not fit/cut) inspection.
//
// Reference: Ghosh, Laramore & McGregor, Nucl. Instrum. Methods Phys.
// Res. A 984 (2020) 164496, Eq. (5):
//   N_n(t) = 0.754*exp(-t/76.37) + 0.1704*exp(-t/702.66) + 0.0754*exp(-t/9920)
// (t in ns), fit to the 100-pulse average neutron-induced pulse, normalized
// so N_n(0) ~= 1 (Fig. 9(a)).
//
// This is a PURE OVERLAY: the curve shape/time constants are fixed from the
// paper, never fit to the event. Only the overall amplitude is scaled so
// the curve's value at t=0 matches the event's own peak (wf.getMaxAmpADC()).
// No accept/reject logic, no chi^2 -- for visual comparison only, following
// from the earlier decision to remove the per-event DecayFit class (a
// chi^2/ndf cut on this same functional form rejected 100% of real events
// due to bumpy real data, not because the underlying shape was wrong).
class NeutronCurveOverlay {
public:
    NeutronCurveOverlay();

    // Draws wf's baseline-subtracted waveform (same style as
    // Waveform::drawWithThreshold()) with the reference neutron decay
    // curve overlaid, anchored at explicitStartSample and scaled to the
    // event's own peak amplitude. Prints to the currently-open multi-page
    // PDF (c and pdfName are managed by the caller, same convention as
    // Waveform::drawWithThreshold()).
    //
    // explicitStartSample: sample index where the curve's t=0 is placed
    // (typically the detected pulse start from
    // Waveform::findNeutronLikeStretch() / countSeparatePulses()).
    // nsPerSample: sample-to-time conversion (same units as elsewhere).
    // extraLabel: optional text drawn in the upper-left corner, same
    // convention as Waveform::drawWithThreshold().
    void draw(const Waveform& wf, int explicitStartSample, double nsPerSample,
              long eventIndex, TCanvas* c, const std::string& pdfName,
              const std::string& extraLabel = "") const;

private:
    // Ghosh, Laramore & McGregor, NIM A 984 (2020) 164496, Eq. (5).
    static constexpr double T_FAST = 76.37;   // ns
    static constexpr double T_MID  = 702.66;  // ns
    static constexpr double T_SLOW = 9920.0;  // ns
    static constexpr double W_FAST = 0.754;
    static constexpr double W_MID  = 0.1704;
    static constexpr double W_SLOW = 0.0754;
    // Weights sum to 0.9998 in the paper (not exactly 1) -- normalize by
    // this so evaluate(0) is exactly 1.0 before external scaling, rather
    // than silently carrying a 0.02% offset.
    static constexpr double W_SUM = W_FAST + W_MID + W_SLOW;

    // Reference curve shape, normalized so evaluate(0) == 1.0. Returns 0
    // for tNs < 0 (curve only defined from the anchor point onward).
    double evaluate(double tNs) const;
};

#endif // NEUTRONCURVEOVERLAY_H