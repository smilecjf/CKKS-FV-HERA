LIBPATH = ../lib
XKCPPATH = $(LIBPATH)/XKCP_AVX2
XKCPHEADERPATH = $(XKCPPATH)/libXKCP.a.headers
XKCPHEADER = $(XKCPHEADERPATH)/KeccakP-1600-times4-SnP.h $(XKCPHEADERPATH)/SIMD256-config.h
SEALPATH = $(LIBPATH)/libSEAL-CKKS-FV
SEALHEADER = $(SEALPATH)/include/SEAL-3.4

CXXOPT = -O3 -I$(XKCPHEADERPATH) -I$(SEALHEADER) -std=c++17
LIBXKCP = -L$(XKCPPATH) -lXKCP
LIBSEAL = -L$(SEALPATH)/lib -lseal-3.4 -lpthread -lm -lz

HERASOURCE = Hera.cpp ShakeAVX2.cpp
HERAHEADER = Hera.h ShakeAVX2.h
ENCODERSOURCE = CKKS_Encoder.cpp FV_Encoder.cpp
ENCODERHEADER = CKKS_Encoder.h FV_Encoder.h


all: sym_enc hera

sym_enc: sym_enc.cpp $(HERASOURCE) $(HERAHEADER) $(ENCODERSOURCE) $(ENCODERHEADER) $(XKCPHEADER) parms.h
	$(CXX) $(CXXOPT) sym_enc.cpp $(HERASOURCE) $(ENCODERSOURCE) $(LIBSEAL) $(LIBXKCP) -o sym_enc -mavx2

hera: test_hera.cpp $(HERASOURCE) $(HERAHEADER) parms.h $(XKCPHEADER)
	$(CXX) $(CXXOPT) test_hera.cpp $(HERASOURCE) $(LIBXKCP) -o test_hera -mavx2


clean:
	rm -f sym_enc
	rm -f test_hera
