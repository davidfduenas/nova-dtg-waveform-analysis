// Utilities.h
#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>

namespace Utilities {

    // Change this to point at a different data location.
    extern const std::string INPUT_DIR;

    // Builds "<macroName>_results/<dataFileBaseName>/" from the macro's
    // own argv[0] and one input filename.
    std::string makeOutputDir(const std::string& argv0, const std::string& inputFileName);

    // Same, for macros taking two input files: builds
    // "<macroName>_results/<baseName1>_vs_<baseName2>/".
    std::string makeOutputDir(const std::string& argv0, const std::string& inputFileName1,
                              const std::string& inputFileName2);

}

#endif // UTILITIES_H
