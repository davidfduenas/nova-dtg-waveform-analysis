// ============================================================================
// Convert.C
// Converts EJ410 raw .txt files into ROOT files (waveforms TTree)
//
// Supports two input formats — auto-detected from the first line of each file:
//
//    WaveDump / CoMPASS:
//     Record Length: 1030
//     BoardID: 31
//     Channel: 0
//     Event Number: 579018
//     Pattern: 0x0000
//     Trigger Time Stamp: 2293227763
//     DC offset (DAC): 0x3333
//     <ADC samples, one per line>
//
//    DAQ threshold (with unixtime, current code version):
//     # Event N  tag=self  trig=self  ch=0  size=1500  cnt=N  ttag=<t>  unixtime=<t>
//     <ADC samples, one per line>
//
//    DT5730 digitizer:
//   - 14-bit ADC, 2Vpp input range
//   - 500 MS/s sampling rate (2 ns per sample)
//   - Timestamp: 31-bit counter + overflow bit, 8 ns increment / 16 ns resolution
//   - Record length is read from each event header
//
// Output ROOT tree:
//   eventNumber           /I   software loop index (resets every DAQ launch)
//   hwCounter             /I   hardware EventCounter (use to check for lost
//                                events between board and file)
//   recordLength          /I
//   triggerTimeStamp      /L   RAW ttag, relative counter, 8 ns/tick
//   unixTime              /D   REAL wall-clock time (seconds since epoch,
//                                with sub-second precision) — use THIS for
//                                any question about when an event actually
//                                happened, or how long a run lasted. Only
//                                present if the .txt file's header includes
//                                "unixtime=" (current DAQ code version); older
//                                files without it will have unixTime=0.
//   waveform[recordLength]/I
//
// Usage:
//   root -l -q 'Convert.C("run.txt")'
// ============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

const Int_t MAX_RECORD_LENGTH = 20000;

// ── format detection ──────────────────────────────────────────────────────────
std::string detectFormat(const std::string& path)
{
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line))
        if (!line.empty()) return (line[0] == '#') ? "new" : "old";
    return "old";
}

// ── OLD FORMAT: read one event ────────────────────────────────────────────────
bool readEventOld(std::ifstream& file, Int_t& recordLength, Int_t& evtNum,
                  Int_t& hwCounter, Long64_t& timestamp, Double_t& unixTime,
                  std::vector<int>& samples)
{
    std::string line;
    while (std::getline(file, line))
        if (line.find("Record Length:") != std::string::npos) {
            recordLength = std::stoi(line.substr(line.find(":") + 1));
            break;
        }
    if (file.eof()) return false;

    std::getline(file, line);  // BoardID
    std::getline(file, line);  // Channel
    std::getline(file, line);  evtNum    = std::stoi(line.substr(line.find(":") + 1));
    hwCounter = evtNum; // no separate hardware counter field in this format
    std::getline(file, line);  // Pattern
    std::getline(file, line);  timestamp = std::stoll(line.substr(line.find(":") + 1));
    unixTime = 0.0; // this format has no wall-clock time field
    std::getline(file, line);  // DC offset

    samples.clear();
    samples.reserve(recordLength);
    for (int i = 0; i < recordLength; i++) {
        if (!std::getline(file, line)) break;
        if (!line.empty()) samples.push_back(std::stoi(line));
    }
    return (int)samples.size() == recordLength;
}

// ── NEW FORMAT: read one event ────────────────────────────────────────────────
bool readEventNew(std::ifstream& file, Int_t& recordLength, Int_t& evtNum,
                  Int_t& hwCounter, Long64_t& timestamp, Double_t& unixTime,
                  std::vector<int>& samples)
{
    std::string line;
    while (std::getline(file, line))
        if (!line.empty() && line[0] == '#') break;
    if (file.eof()) return false;

    evtNum = 0; hwCounter = 0; recordLength = 0; timestamp = 0; unixTime = 0.0;

    {
        std::istringstream ss(line);
        std::string tok;
        ss >> tok >> tok >> evtNum;  // '#' 'Event' N
    }

    std::istringstream ss(line);
    std::string tok;
    while (ss >> tok) {
        if (tok.substr(0, 5) == "size=")     recordLength = std::stoi(tok.substr(5));
        if (tok.substr(0, 4) == "cnt=")      hwCounter    = std::stoi(tok.substr(4));
        if (tok.substr(0, 5) == "ttag=")     timestamp    = std::stoll(tok.substr(5));
        if (tok.substr(0, 9) == "unixtime=") unixTime     = std::stod(tok.substr(9));
    }
    if (recordLength <= 0) return false;

    samples.clear();
    samples.reserve(recordLength);
    for (int i = 0; i < recordLength; i++) {
        if (!std::getline(file, line)) break;
        if (!line.empty() && line[0] != '#') samples.push_back(std::stoi(line));
        else --i;
    }
    return (int)samples.size() == recordLength;
}

// ── unified reader ────────────────────────────────────────────────────────────
bool readEvent(std::ifstream& file, const std::string& fmt,
               Int_t& recordLength, Int_t& evtNum, Int_t& hwCounter,
               Long64_t& timestamp, Double_t& unixTime, std::vector<int>& samples)
{
    return (fmt == "new")
        ? readEventNew(file, recordLength, evtNum, hwCounter, timestamp, unixTime, samples)
        : readEventOld(file, recordLength, evtNum, hwCounter, timestamp, unixTime, samples);
}

// ── convert one file ──────────────────────────────────────────────────────────
void convertOneFile(const std::string& path)
{
    std::string fmt = detectFormat(path);

    std::string rootPath = path;
    size_t dot = rootPath.rfind(".txt");
    if (dot != std::string::npos) rootPath.replace(dot, 4, ".root");
    else rootPath += ".root";

    std::ifstream infile(path);
    if (!infile.is_open()) {
        std::cerr << "WARNING: cannot open " << path << std::endl;
        return;
    }

    gSystem->Unlink(rootPath.c_str());
    TFile fOut(rootPath.c_str(), "RECREATE");

    Int_t    evtNum, hwCounter, recordLength;
    Int_t    waveform[MAX_RECORD_LENGTH];
    Long64_t timestamp;
    Double_t unixTime;

    TTree tree("waveforms", "EJ410 detector data");
    tree.Branch("eventNumber",      &evtNum,       "eventNumber/I");
    tree.Branch("hwCounter",        &hwCounter,    "hwCounter/I");
    tree.Branch("recordLength",     &recordLength, "recordLength/I");
    tree.Branch("triggerTimeStamp", &timestamp,    "triggerTimeStamp/L");
    tree.Branch("unixTime",         &unixTime,     "unixTime/D");
    tree.Branch("waveform",          waveform,     "waveform[recordLength]/I");

    Long64_t n = 0;
    Int_t lastHwCounter = -1;
    Long64_t gapCount = 0;
    std::vector<int> samples;
    while (readEvent(infile, fmt, recordLength, evtNum, hwCounter, timestamp, unixTime, samples)) {
        for (int i = 0; i < recordLength; i++) waveform[i] = samples[i];

        if (lastHwCounter >= 0 && hwCounter != lastHwCounter + 1) {
            gapCount++;
            std::cerr << "  [warn] hwCounter discontinuity: " << lastHwCounter
                      << " -> " << hwCounter << " (row " << n
                      << ", unixTime=" << std::fixed << unixTime << ")\n";
        }
        lastHwCounter = hwCounter;

        tree.Fill(); n++;
        if (n % 10000 == 0)
            std::cout << "  " << path << ": " << n << " events\r" << std::flush;
    }
    tree.Write("", TObject::kOverwrite);
    fOut.Close(); infile.close();

    std::cout << "\n  " << path << " -> " << rootPath
              << "  (" << n << " events, recLen=" << recordLength
              << ", fmt=" << fmt << ", counter gaps=" << gapCount << ")" << std::endl;
}

// ── entry point ───────────────────────────────────────────────────────────────
void Convert(const char* files = "")
{
    // >>> SET YOUR DATA DIRECTORY HERE <
    std::string dataPath = "/Users/david/DTGdata/EJ410data/testruns/";

    std::vector<std::string> fileList;
    std::istringstream ss(files);
    std::string tok;
    while (ss >> tok) fileList.push_back(dataPath + tok);

    if (fileList.empty()) {
        std::cerr << "Usage: Convert(\"file1.txt [file2.txt ...]\")\n";
        return;
    }

    std::cout << "\n=== Converting files ============================\n" << std::endl;
    for (const auto& f : fileList)
        convertOneFile(f);
    std::cout << "\n=== Done ========================================\n" << std::endl;
}