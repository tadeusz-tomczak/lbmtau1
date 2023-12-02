#
# Uncomment one of the options below.
#
OPTIONS = -DKERNEL_OPT_MANUAL
#OPTIONS = -DKERNEL_VECTORIZED_LOOP
#OPTIONS = -DKERNEL_REFERENCE
#OPTIONS = -DKERNEL_COPY_ONLY

OPTIONS += -DUSE_LARGE_PAGES 

INCLUDES = -I.
LIBS     = -lstdc++ -lm

CLANG = clang $(OPTIONS) $(EXTRA_OPTIONS) -march=native -Ofast -zopt -freciprocal-math -mrecip -finline-aggressive -faggressive-loop-transform  -ffast-math -floop-reroll -floop-splitting -floop-transform  -flto  -fdebug-info-for-profiling -g  -fopenmp --std=c++17 $(CXXFLAGS) -Wall -pthread $(INCLUDES)

clang: CPP = $(CLANG)
clang: main

main: main.cpp LBMTau1.hpp
	$(CPP) main.cpp -o $@ $(LIBS)


.PHONY: clean
clean:
	rm -rf main

