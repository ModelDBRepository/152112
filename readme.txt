readme.txt

Simulation code accompanying the manuscript:
"A model of the medial superior olive explains spatiotemporal features
of local field potentials"
By JH Goldwyn, M McLaughlin, E Verschooten, PX Joris, J Rinzel
Manuscript submitted to J Neuroscience 1/14/2014

Matlab (R2012b) simulation code written by JH Goldwyn and posted to
ModelDB on 1/14/2014

This code makes use of SUNDIALS (Suite of Nonlinear and Differential
Algebraic Equation Solvers) and its interface to Matlab (sundialsTB).

These can be downloaded at the website:

http://computation.llnl.gov/casc/sundials/main.html

Documentation and installation instructions for SUNDIALS and
sundialsTB are also available at that address.

Contents:
makeFig.m: An m-file that reproduces Figures 4, 8, and 11 from the
manuscript.

runNeurophonic.m: An m-file that produces 3D (surf command)
spatial-temporal profiles of membrane potential and extracellular
potential.  Parameter values in this file that can be modified by the
user include synaptic strength (gE, gI); frequency of synaptic events
(synFreq); dendrite receiving excitatory synaptic events (stimType);
length of simulation (tEnd); etc.

MSO_dae.m: A function file that defines and solves the system of
equations that model the membrane potential of a MSO neuron and
extracellular potential in the surrounding "virtual cylinder" of
extracellular space.  See manuscript for details. This function is
called by makeFig.m and runNeurophonic.m
