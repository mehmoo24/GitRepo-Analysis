
#################################################
# Makefile automatically generated by:
# ../scripts/GenerateMakefile.py --dir .
#################################################

CXX = g++ 
CXXFLAGS = -g -Wall -fPIC 
ROOTFLAGS = `root-config --cflags --glibs --libs`
#INCLUDE += -I$(NUKECCSRC_INCLUDE_PATH) -I$(NUKECCSRC_ANA) -I$(PLOTUTILSROOT/../include)/ -I$(UNFOLDUTILSROOT)/
#LDLIBS += -L$(PLOTUTILSROOT)/ -lplotutils -L$(UNFOLDUTILSROOT)/$(CMTCONFIG) -lUnfoldUtils
INCLUDE += -I$(NUKECCSRC_INCLUDE_PATH) -I$(NUKECCSRC_ANA) -I$(PLOTUTILSROOT)/../include/
LDLIBS += -L$(PLOTUTILSROOT) -lMAT -lMAT-MINERvA -lGenVector 

LDLIBS += -L$(NUKECCSRC_ANA)/src -lNukeCC -lUnfoldUtils

LDLIBS += -L$(GENIEXSECEXTRACTROOT) -lGENIEXSecExtract

TARGETS   = $(wildcard *.o)
BINARIES  = $(TARGETS:.o=)
all : $(patsubst %.cxx, %.o, $(wildcard *.cxx))

#--- if making individual targets ---#                                                                    
%.o : %.cxx $(DEPS)
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) -o $@ $(LDLIBS) -c $< #compile                         
	$(CXX) $(INCLUDE) $(CXXFLAGS) $(ROOTFLAGS) $(LDLIBS) -o $* $@        #link                        
.PHONY: clean

clean:
	rm -f $(BINARIES) $(TARGETS)
