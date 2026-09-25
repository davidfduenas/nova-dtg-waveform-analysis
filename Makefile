ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

SRC_DIR    := src
LIB_DIR    := lib
MACRO_DIR  := macros
BUILD_DIR  := build

LIB_TARGET := $(LIB_DIR)/libWaveform.so

# Automatically find every .cpp file in macros/, and build a matching
# executable (same name, no extension) in build/ for each one.
MACRO_SOURCES := $(wildcard $(MACRO_DIR)/*.cpp)
MACRO_TARGETS := $(patsubst $(MACRO_DIR)/%.cpp,$(BUILD_DIR)/%,$(MACRO_SOURCES))

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LIB_INSTALL_NAME := -install_name @rpath/libWaveform.so
    EXE_RPATH        := -Wl,-rpath,@executable_path/../lib
else
    LIB_INSTALL_NAME :=
    EXE_RPATH        := -Wl,-rpath,'$$ORIGIN/../lib'
endif

.PHONY: all clean

all: $(LIB_TARGET) $(MACRO_TARGETS)

$(LIB_TARGET): $(SRC_DIR)/Waveform.cpp $(SRC_DIR)/Waveform.h $(SRC_DIR)/PSD.cpp $(SRC_DIR)/PSD.h $(SRC_DIR)/NeutronCurveOverlay.cpp $(SRC_DIR)/NeutronCurveOverlay.h
	mkdir -p $(LIB_DIR)
	g++ -O2 -shared -fPIC $(ROOTCFLAGS) $(SRC_DIR)/Waveform.cpp $(SRC_DIR)/PSD.cpp $(SRC_DIR)/NeutronCurveOverlay.cpp -o $(LIB_TARGET) $(ROOTLIBS) $(LIB_INSTALL_NAME)

# Pattern rule: "to build build/<name>, from macros/<name>.cpp"
$(BUILD_DIR)/%: $(MACRO_DIR)/%.cpp $(LIB_TARGET)
	mkdir -p $(BUILD_DIR)
	g++ -O2 $(ROOTCFLAGS) -I$(SRC_DIR) $< \
	    -L$(LIB_DIR) -lWaveform $(ROOTLIBS) $(EXE_RPATH) \
	    -o $@

clean:
	rm -f $(LIB_TARGET)
	rm -rf $(BUILD_DIR)