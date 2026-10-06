ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

SRC_DIR    := src
LIB_DIR    := lib
MACRO_DIR  := macros

LIB_TARGET := $(LIB_DIR)/libWaveform.so

# Every .cpp in src/ goes into the shared library.
LIB_SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
LIB_HEADERS := $(wildcard $(SRC_DIR)/*.h)

# Automatically find every .cpp file in macros/, and build a matching
# executable (same name, no extension) for each one.
MACRO_SOURCES := $(wildcard $(MACRO_DIR)/*.cpp) $(wildcard $(MACRO_DIR)/*/*.cpp)
MACRO_TARGETS := $(MACRO_SOURCES:.cpp=)

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LIB_INSTALL_NAME := -install_name @rpath/libWaveform.so
    EXE_RPATH        := -Wl,-rpath,@executable_path/../../lib
else
    LIB_INSTALL_NAME :=
    EXE_RPATH        := -Wl,-rpath,'$$ORIGIN/../../lib'
endif

.PHONY: all clean update-gitignore

all: $(LIB_TARGET) $(MACRO_TARGETS) update-gitignore

# Rewrites just the marked section of .gitignore with the CURRENT exact
# list of compiled macro binaries (MACRO_TARGETS). Matches filenames
# exactly, never a wildcard, so it can't accidentally ignore some other
# output type later. Runs automatically after every `make all`, so a
# newly added macro's binary is picked up with no manual edit needed.
update-gitignore:
	@awk -v list="$(MACRO_TARGETS)" ' \
	    /# BEGIN AUTO-GENERATED MACRO BINARIES/ { \
	        print; \
	        n = split(list, arr, " "); \
	        for (i = 1; i <= n; i++) print arr[i]; \
	        inblock = 1; \
	        next \
	    } \
	    /# END AUTO-GENERATED MACRO BINARIES/ { inblock = 0 } \
	    !inblock { print } \
	' .gitignore > .gitignore.tmp && mv .gitignore.tmp .gitignore

$(LIB_TARGET): $(LIB_SOURCES) $(LIB_HEADERS)
	mkdir -p $(LIB_DIR)
	g++ -O2 -shared -fPIC $(ROOTCFLAGS) $(LIB_SOURCES) -o $(LIB_TARGET) $(ROOTLIBS) $(LIB_INSTALL_NAME)

# Pattern rule: "to build any file with no extension, from a matching
# .cpp file of the same name"
$(MACRO_DIR)/%: $(MACRO_DIR)/%.cpp $(LIB_TARGET)
	g++ -O2 $(ROOTCFLAGS) -I$(SRC_DIR) $< \
	    -L$(LIB_DIR) -lWaveform $(ROOTLIBS) $(EXE_RPATH) \
	    -o $@

clean:
	rm -f $(LIB_TARGET) $(MACRO_TARGETS)