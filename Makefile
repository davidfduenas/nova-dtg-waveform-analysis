ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS   := $(shell root-config --libs)

SRC_DIR    := src
LIB_DIR    := lib
MACRO_DIR  := macros

LIB_TARGET := $(LIB_DIR)/libWaveform.so

# Automatically find every .cpp file in macros/, and build a matching
# executable (same name, no extension) for each one.
MACRO_SOURCES := $(wildcard $(MACRO_DIR)/*.cpp)
MACRO_TARGETS := $(MACRO_SOURCES:.cpp=)

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

$(LIB_TARGET): $(SRC_DIR)/Waveform.cpp $(SRC_DIR)/Waveform.h
	mkdir -p $(LIB_DIR)
	g++ -O2 -shared -fPIC $(ROOTCFLAGS) $(SRC_DIR)/Waveform.cpp -o $(LIB_TARGET) $(ROOTLIBS) $(LIB_INSTALL_NAME)

# Pattern rule: "to build any file with no extension, from a matching
# .cpp file of the same name" -- this one rule replaces having to write
# a separate rule for every single program.
$(MACRO_DIR)/%: $(MACRO_DIR)/%.cpp $(LIB_TARGET)
	g++ -O2 $(ROOTCFLAGS) -I$(SRC_DIR) $< \
	    -L$(LIB_DIR) -lWaveform $(ROOTLIBS) $(EXE_RPATH) \
	    -o $@

clean:
	rm -f $(LIB_TARGET) $(MACRO_TARGETS)