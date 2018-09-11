ALLCPP2 = $(wildcard */*.cpp) 
ALLCPP = $(ALLCPP2:boost_tests/%=) 
ALLH = $(wildcard */*.h) 
VER = -DVERSION=$(shell cat VERSION)
BOOST_PO = -lboost_program_options
BOOST_TH = -lboost_thread
BOOST_S = -lboost_system
PTHREAD= -lpthread
NOBAMWAR= -std=c++0x -Wno-unknown-pragmas
Z = -lz
M = -lm
BOOST_PATH = 
STATIC =

all:peakranger  
peakranger    : $(ALLCPP) $(ALLH) 
	mkdir -p bin
	g++ -I./ $(BOOST_PATH) $(STATIC) $(NOBAMWAR) -Wall -pedantic  -O3  $(ALLCPP)  -o bin/$@ $(VER) $(BOOST_S) $(BOOST_PO) $(BOOST_FS) $(Z) $(M) $(BOOST_TH)  $(PTHREAD) -std=c++0x
	@echo "PeakRanger compilation complete."
%.o : %.cpp
	g++ -Wall -pedantic $(NOBAMWAR) $(LOG) $(VER) -O3  -c $< -o $@
