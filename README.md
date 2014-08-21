# SpikeSorter
This is a package for sorting so-called spikes extracted from electrophysiological data recorded extra-cellularly. The main functionality is to assign combinations of spike templates to 
high pass filtered data using the viterbi algorithms. The code is based on the following paper:
Herbst, J. A., Gammeter, S., Ferrero, D. & Hahnloser, R. H. R. Spike sorting with hidden Markov models. Journal of Neuroscience Methods 174, 126–134 (2008).

##Introduction

Spike sorting is the art of decomposing a signal recorded from the extra-cellular medium into contributions from individual cells. These contributions vary in shape depending on a cell's membrane potential dynamics, as well as the location of that cell in relation to the recording probe. Needless to say, this is a serverly under constrained problem, and various approximations need to be made in order to settle on a plausible solution.

###Hidden markov model

In Hidden Markov model sorting, the part of the extra-cellular signal that can be attributed to a single cells is described as a set of state transitions with associated gaussian emission probablities.

[![Build Status](https://travis-ci.org/roger.herikstad@gmail.com/SpikeSorter.jl.png)](https://travis-ci.org/roger.herikstad@gmail.com/SpikeSorter.jl)
