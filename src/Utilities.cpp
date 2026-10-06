// Utilities.cpp
#include "Utilities.h"

namespace Utilities {

const std::string INPUT_DIR = "/Users/david/DTGAnalysis/data/collruns/";

static std::string baseName(const std::string& path)
{
    std::string name = path;
    size_t slash = name.find_last_of('/');
    if (slash != std::string::npos) name = name.substr(slash + 1);
    size_t dot = name.find_last_of('.');
    if (dot != std::string::npos) name = name.substr(0, dot);
    return name;
}

static std::string macroName(const std::string& argv0)
{
    std::string name = argv0;
    size_t slash = name.find_last_of('/');
    if (slash != std::string::npos) name = name.substr(slash + 1);
    return name;
}

std::string makeOutputDir(const std::string& argv0, const std::string& inputFileName)
{
    return macroName(argv0) + "_results/" + baseName(inputFileName) + "/";
}

std::string makeOutputDir(const std::string& argv0, const std::string& inputFileName1,
                          const std::string& inputFileName2)
{
    return macroName(argv0) + "_results/" + baseName(inputFileName1) + "_vs_" + baseName(inputFileName2) + "/";
}

}
