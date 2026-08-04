// ============================================================================
// PulseFinder.C
// Finds ALL pulses in each waveform — not just the largest one
//
// For each event, scans the waveform and identifies every region where
// the signal dips below the per-event sigma-cut baseline. Each dip is
// recorded as a separate pulse with its own height and area.
//
// This allows:
//   - Seeing the full pulse height spectrum including small pulses
//   - Separating gammas/x-rays (small) from neutrons (large)
//   - Comparing multiple runs directly
//
// Detection:
//   A pulse is a contiguous run of samples dipping more than START_K·sigma
//   below baseline. Must be at least MIN_WIDTH samples wide.
//
// Output directories are derived from input filenames automatically.
// No filenames or paths are hardcoded.
//
// Usage:
//   // Single file:
//   root -l -q 'PulseFinder.C("run.root")'
//
//   // Two files (overlay):
//   root -l -q 'PulseFinder.C("source.root", "background.root")'
//
//   // Three files (three-way overlay):
//   root -l -q 'PulseFinder.C("background.root", "source.root", "container.root")'
// ============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TSystem.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// ── digitizer constants ───────────────────────────────────────────────────────
const double NS_PER_SAMPLE     = 2.0;
const double NS_PER_TICK       = 8.0;
const int    MAX_RECORD_LENGTH = 20000;

// ── sigma-cut baseline dials ──────────────────────────────────────────────────
const double SIGMA_CUT_A        = 3.0;
const int    SIGMA_CUT_MAX_ITER = 10;
const double SIGMA_CUT_TOL      = 0.01;

// ── pulse-detection dials ─────────────────────────────────────────────────────
const double START_K   = 5.0;  // threshold = START_K · sigma below baseline
const int    MIN_WIDTH = 1;    // minimum pulse width in samples

// ── derive output directory name from filename ────────────────────────────────
std::string outDir(const std::string& fileName)
{
    size_t slash = fileName.rfind('/');
    std::string base = (slash == std::string::npos) ? fileName : fileName.substr(slash + 1);
    size_t dot = base.rfind(".root");
    if (dot != std::string::npos) base = base.substr(0, dot);
    return "PulseFinder_" + base;
}

// ── create output subdirectory tree ──────────────────────────────────────────
void makeOutDirs(const std::string& dir)
{
    gSystem->mkdir(dir.c_str(), true);
    gSystem->mkdir((dir + "/tagged_waveforms").c_str(), true);
    gSystem->mkdir((dir + "/spectra").c_str(),          true);
    gSystem->mkdir((dir + "/multiplicity").c_str(),     true);
}

// ── style ─────────────────────────────────────────────────────────────────────
void setStyle() {
    gStyle->SetOptStat(0); gStyle->SetOptTitle(0);
    gStyle->SetPadGridX(1); gStyle->SetPadGridY(1);
    gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetFrameLineWidth(2); gStyle->SetHistLineWidth(2);
    gStyle->SetLabelSize(0.045, "XY"); gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(1.1, "Y");
}

// ── sigma-cut baseline ────────────────────────────────────────────────────────
double computeBaselineSigmaCut(const Int_t* wf, int n,
                               double A = SIGMA_CUT_A, int maxIter = SIGMA_CUT_MAX_ITER,
                               double tol = SIGMA_CUT_TOL, double* sigmaOut = nullptr)
{
    double sum = 0, sum2 = 0;
    for (int i = 0; i < n; i++) { sum += wf[i]; sum2 += (double)wf[i]*wf[i]; }
    double mean = sum/n, var = (sum2 - n*mean*mean)/(n-1);
    double sigma = (var > 0) ? std::sqrt(var) : 0;

    for (int iter = 0; iter < maxIter; iter++) {
        if (sigma <= 0) break;
        double lo = mean - A*sigma, hi = mean + A*sigma;
        sum = 0; sum2 = 0; long m = 0;
        for (int i = 0; i < n; i++)
            if (wf[i] >= lo && wf[i] <= hi) { sum += wf[i]; sum2 += (double)wf[i]*wf[i]; m++; }
        if (m < 2) break;
        double newMean = sum/m;
        var = (sum2 - m*newMean*newMean)/(m-1);
        double newSigma = (var > 0) ? std::sqrt(var) : 0;
        bool converged = std::fabs(newMean - mean) < tol;
        mean = newMean; sigma = newSigma;
        if (converged) break;
    }
    if (sigmaOut) *sigmaOut = sigma;
    return mean;
}

// ── pulse structure ───────────────────────────────────────────────────────────
struct Pulse {
    double height, area;
    int start, end, peak, width;
};

// ── find all pulses in a waveform ─────────────────────────────────────────────
std::vector<Pulse> findPulses(const Int_t* waveform, int recordLength,
                              double baseline, double sigma,
                              double startK = START_K, int minWidth = MIN_WIDTH,
                              int* nRejected = nullptr) {
    std::vector<Pulse> pulses;
    double thresh = startK * sigma;
    int rejected = 0;
    bool inPulse = false;
    Pulse current;

    auto closePulse = [&](int endIdx) {
        current.end   = endIdx;
        current.width = current.end - current.start + 1;
        if (current.width >= minWidth) pulses.push_back(current);
        else rejected++;
        inPulse = false;
    };

    for (int i = 0; i < recordLength; i++) {
        double dev = baseline - waveform[i];
        if (dev > thresh) {
            if (!inPulse) {
                inPulse = true;
                current.start = i; current.height = dev; current.area = dev; current.peak = i;
            } else {
                current.area += dev;
                if (dev > current.height) { current.height = dev; current.peak = i; }
            }
        } else if (inPulse) closePulse(i - 1);
    }
    if (inPulse) closePulse(recordLength - 1);
    if (nRejected) *nRejected = rejected;
    return pulses;
}

// ── run duration from timestamps ──────────────────────────────────────────────
double getDuration(TTree* tree) {
    return (tree->GetMaximum("triggerTimeStamp") -
            tree->GetMinimum("triggerTimeStamp")) * NS_PER_TICK * 1e-9;
}

// ── plot tagged example waveforms ────────────────────────────────────────────
void plotTaggedWaveforms(const char* fileName, const char* label,
                         const char* prefix, int nPlots = 6, int minPulses = 1) {
    TFile* f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fileName << std::endl; return; }
    TTree* tree = (TTree*)f->Get("waveforms");
    if (!tree) { f->Close(); return; }

    Int_t waveform[MAX_RECORD_LENGTH], recordLength;
    tree->SetBranchAddress("waveform",     waveform);
    tree->SetBranchAddress("recordLength", &recordLength);

    Long64_t nEntries = tree->GetEntries();
    int made = 0;

    for (Long64_t i = 0; i < nEntries && made < nPlots; i++) {
        tree->GetEntry(i);
        double sigma;
        double baseline = computeBaselineSigmaCut(waveform, recordLength,
                                                  SIGMA_CUT_A, SIGMA_CUT_MAX_ITER,
                                                  SIGMA_CUT_TOL, &sigma);
        std::vector<Pulse> pulses = findPulses(waveform, recordLength, baseline, sigma);
        if ((int)pulses.size() < minPulses) continue;

        std::cout << "  [" << label << " evt " << i << "] sigma=" << Form("%.2f",sigma)
                  << " thresh=" << Form("%.1f",START_K*sigma) << " ADC, "
                  << pulses.size() << " pulses" << std::endl;
        for (size_t pi = 0; pi < pulses.size(); pi++) {
            const Pulse& p = pulses[pi];
            std::cout << "    pulse " << pi << ": height=" << Form("%.1f",p.height)
                      << "  peak@" << p.peak << " (" << p.peak*NS_PER_SAMPLE << "ns)"
                      << "  width=" << p.width << std::endl;
        }

        double windowNs = recordLength * NS_PER_SAMPLE;
        TGraph* gr = new TGraph(recordLength);
        for (int s = 0; s < recordLength; s++)
            gr->SetPoint(s, s * NS_PER_SAMPLE, waveform[s] - baseline);

        TCanvas* c = new TCanvas(Form("c%s%d", prefix, made),
                                  Form("%s tagged %d", label, made+1), 1200, 600);
        c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
        gr->GetXaxis()->SetTitle("Time [ns]");
        gr->GetYaxis()->SetTitle("ADC Counts");
        gr->SetLineColor(kBlue+1); gr->SetLineWidth(2);
        gr->Draw("AL");

        TLine* lBase = new TLine(0, 0, windowNs, 0);
        lBase->SetLineColor(kRed); lBase->SetLineWidth(2); lBase->SetLineStyle(2); lBase->Draw();

        for (const Pulse& p : pulses) {
            TMarker* mk = new TMarker(p.peak * NS_PER_SAMPLE, -p.height, 32);
            mk->SetMarkerColor(kRed); mk->SetMarkerSize(1.0); mk->Draw();
            TLine* lS = new TLine(p.start * NS_PER_SAMPLE, 0, p.start * NS_PER_SAMPLE, -p.height);
            lS->SetLineColor(kGray+1); lS->SetLineStyle(3); lS->Draw();
            TLine* lE = new TLine((p.end+1) * NS_PER_SAMPLE, 0, (p.end+1) * NS_PER_SAMPLE, -p.height);
            lE->SetLineColor(kGray+1); lE->SetLineStyle(3); lE->Draw();
        }

        TLatex* title = new TLatex(0.5, 0.95,
            Form("EJ-410 %s (event %lld, %d pulses)", label, i, (int)pulses.size()));
        title->SetNDC(); title->SetTextAlign(22); title->SetTextSize(0.045); title->Draw();

        c->SaveAs(Form("%s_%02d.pdf", prefix, made+1));
        std::cout << "Saved " << prefix << "_"
                  << std::setw(2) << std::setfill('0') << made+1 << ".pdf"
                  << " (event " << i << ", " << pulses.size() << " pulses)" << std::endl;
        made++;
    }
    if (made == 0)
        std::cout << "  (no events with >= " << minPulses << " pulses in " << fileName << ")" << std::endl;
    f->Close();
}

// ── process one file: fill histograms ────────────────────────────────────────
Long64_t processFile(const char* fileName, TH1D* hHeight, TH1D* hArea,
                     TH1D* hMult, double& durationOut, Long64_t* nEventsOut = nullptr) {
    TFile* f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fileName << std::endl; durationOut = 0; return 0; }
    TTree* tree = (TTree*)f->Get("waveforms");
    if (!tree) { f->Close(); durationOut = 0; return 0; }

    Int_t waveform[MAX_RECORD_LENGTH], recordLength;
    tree->SetBranchAddress("waveform",     waveform);
    tree->SetBranchAddress("recordLength", &recordLength);
    durationOut = getDuration(tree);

    Long64_t nEntries = tree->GetEntries(), totalPulses = 0, totalRejected = 0;
    double sumSigma = 0;

    for (Long64_t i = 0; i < nEntries; i++) {
        tree->GetEntry(i);
        double sigma;
        double baseline = computeBaselineSigmaCut(waveform, recordLength,
                                                  SIGMA_CUT_A, SIGMA_CUT_MAX_ITER,
                                                  SIGMA_CUT_TOL, &sigma);
        sumSigma += sigma;
        int nRej = 0;
        std::vector<Pulse> pulses = findPulses(waveform, recordLength, baseline, sigma,
                                               START_K, MIN_WIDTH, &nRej);
        totalRejected += nRej;
        hMult->Fill(pulses.size());
        for (const Pulse& p : pulses) { hHeight->Fill(p.height); hArea->Fill(p.area); totalPulses++; }
        if (i % 10000 == 0) std::cout << "  " << fileName << ": " << i << " events\r" << std::flush;
    }

    double meanSigma = nEntries > 0 ? sumSigma/nEntries : 0;
    std::cout << "\n  " << fileName << ": " << nEntries << " events, "
              << totalPulses << " pulses, duration=" << durationOut << "s"
              << "\n    mean noise sigma=" << Form("%.3f", meanSigma)
              << " ADC  threshold=" << Form("%.2f", START_K*meanSigma) << " ADC"
              << "\n    rejected (width<" << MIN_WIDTH << "): " << totalRejected << std::endl;

    f->Close();
    if (nEventsOut) *nEventsOut = nEntries;
    return totalPulses;
}

// ── draw overlay of up to three histograms ────────────────────────────────────
void drawOverlay(std::vector<TH1D*> hists, std::vector<std::string> labels,
                 const char* xTitle, const char* yTitle,
                 const char* titleText, const char* outFile, bool logy = true) {
    std::vector<int> colors = {kBlack, kRed, kBlue+1, kGreen+2};
    TCanvas* c = new TCanvas(Form("cov_%s", outFile), titleText, 1000, 700);
    c->SetLeftMargin(0.12); c->SetBottomMargin(0.12);
    if (logy) c->SetLogy();

    TLegend* leg = new TLegend(0.55, 0.72, 0.88, 0.88);
    leg->SetBorderSize(1);

    for (size_t i = 0; i < hists.size(); i++) {
        hists[i]->SetLineColor(colors[i % colors.size()]);
        hists[i]->SetLineWidth(2);
        hists[i]->GetXaxis()->SetTitle(xTitle);
        hists[i]->GetYaxis()->SetTitle(yTitle);
        hists[i]->Draw(i == 0 ? "HIST" : "HIST SAME");
        leg->AddEntry(hists[i], labels[i].c_str(), "l");
    }
    leg->Draw();
    TLatex* t = new TLatex(0.5, 0.95, titleText);
    t->SetNDC(); t->SetTextAlign(22); t->SetTextSize(0.045); t->Draw();
    c->SaveAs(outFile);
    std::cout << "Saved " << outFile << std::endl;
}

// ── analyze one file ──────────────────────────────────────────────────────────
void analyzeOne(const char* fileName)
{
    std::string dir = outDir(fileName);
    makeOutDirs(dir);
    std::string label = outDir(fileName).substr(12); // strip "PulseFinder_"

    std::cout << "\n=== PulseFinder: " << fileName << " -> " << dir << "/ ===\n"
              << "  Detection threshold = " << START_K << " sigma,  min width = "
              << MIN_WIDTH << " samples\n" << std::endl;

    // tagged waveforms
    plotTaggedWaveforms(fileName, label.c_str(),
                        (dir + "/tagged_waveforms/wf").c_str(), 6, 1);

    // histograms
    TH1D* hHeight = new TH1D(Form("hH_%s", label.c_str()), "", 200, 0, 600);
    TH1D* hArea   = new TH1D(Form("hA_%s", label.c_str()), "", 200, 0, 5000);
    TH1D* hMult   = new TH1D(Form("hM_%s", label.c_str()), "", 20, 0, 20);

    double dur; Long64_t nEv;
    processFile(fileName, hHeight, hArea, hMult, dur, &nEv);
    if (dur <= 0) { std::cerr << "ERROR: zero duration for " << fileName << std::endl; return; }

    hHeight->Scale(1.0/dur); hArea->Scale(1.0/dur); hMult->Scale(1.0/dur);

    drawOverlay({hHeight}, {label},
                "Pulse Height [ADC counts]", "Pulses / s / bin",
                Form("EJ-410 Pulse Height Rate - %s", label.c_str()),
                (dir + "/spectra/pulse_heights_rate.pdf").c_str());

    drawOverlay({hArea}, {label},
                "Pulse Area [ADC counts #times samples]", "Pulses / s / bin",
                Form("EJ-410 Pulse Area Rate - %s", label.c_str()),
                (dir + "/spectra/pulse_areas_rate.pdf").c_str());

    drawOverlay({hMult}, {label},
                "Number of pulses per event", "Events / s",
                Form("EJ-410 Pulse Multiplicity - %s", label.c_str()),
                (dir + "/multiplicity/multiplicity_rate.pdf").c_str(), false);
}

// ── analyze two or three files with overlays ──────────────────────────────────
void analyzeMulti(std::vector<std::string> fileNames)
{
    int N = fileNames.size();

    // analyze each individually
    for (const auto& f : fileNames) analyzeOne(f.c_str());

    // build shared output dir name
    std::string sharedDir = "PulseFinder";
    for (const auto& f : fileNames) sharedDir += "_" + outDir(f).substr(12);
    gSystem->mkdir(sharedDir.c_str(), true);
    gSystem->mkdir((sharedDir + "/spectra").c_str(), true);
    gSystem->mkdir((sharedDir + "/multiplicity").c_str(), true);

    std::vector<std::string> labels;
    for (const auto& f : fileNames) labels.push_back(outDir(f).substr(12));

    std::cout << "\n=== Multi-file overlays: " << sharedDir << "/ ===\n" << std::endl;

    std::vector<TH1D*> hHeights, hAreas, hMults;
    std::vector<double> durs;

    for (int i = 0; i < N; i++) {
        TH1D* hH = new TH1D(Form("hHm_%d", i), "", 200, 0, 600);
        TH1D* hA = new TH1D(Form("hAm_%d", i), "", 200, 0, 5000);
        TH1D* hM = new TH1D(Form("hMm_%d", i), "", 20, 0, 20);
        double dur; Long64_t nEv;
        processFile(fileNames[i].c_str(), hH, hA, hM, dur, &nEv);
        if (dur > 0) { hH->Scale(1.0/dur); hA->Scale(1.0/dur); hM->Scale(1.0/dur); }
        hHeights.push_back(hH); hAreas.push_back(hA); hMults.push_back(hM);
        durs.push_back(dur);
    }

    drawOverlay(hHeights, labels,
                "Pulse Height [ADC counts]", "Pulses / s / bin",
                "EJ-410 Pulse Height Rate",
                (sharedDir + "/spectra/pulse_heights_rate.pdf").c_str());

    drawOverlay(hAreas, labels,
                "Pulse Area [ADC counts #times samples]", "Pulses / s / bin",
                "EJ-410 Pulse Area Rate",
                (sharedDir + "/spectra/pulse_areas_rate.pdf").c_str());

    drawOverlay(hMults, labels,
                "Number of pulses per event", "Events / s",
                "EJ-410 Pulse Multiplicity Rate",
                (sharedDir + "/multiplicity/multiplicity_rate.pdf").c_str(), false);
}

// ── entry point ───────────────────────────────────────────────────────────────
void PulseFinder(const char* fileA = "",
                 const char* fileB = "",
                 const char* fileC = "")
{
    if (std::string(fileA).empty()) {
        std::cerr << "Usage:\n"
                  << "  PulseFinder(\"file.root\")                          // single file\n"
                  << "  PulseFinder(\"fileA.root\", \"fileB.root\")           // two files\n"
                  << "  PulseFinder(\"fileA.root\", \"fileB.root\", \"fileC.root\") // three files\n";
        return;
    }

    setStyle();

    bool hasB = !std::string(fileB).empty();
    bool hasC = !std::string(fileC).empty();

    if (!hasB && !hasC) {
        analyzeOne(fileA);
    } else {
        std::vector<std::string> files = {fileA};
        if (hasB) files.push_back(fileB);
        if (hasC) files.push_back(fileC);
        analyzeMulti(files);
    }

    std::cout << "\n=== Done ========================================\n" << std::endl;
}
