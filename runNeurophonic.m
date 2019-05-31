%%% THIS CODE CALLS MSO_dae.m %%%

% Computes Vm of MSO neuron model and Ve in extracellular "virtual cylinder"
% Accompanies the manuscript [submitted to J. Neuroscience]:
% "A model of the medial superior olive explains spatiotemporal features of local field potentials"
% JH Goldwyn, M Mc Laughlin, E Verschooten, PX Joris, J Rinzel
%
% MSO model neuron was introduced in:
% "Control of submillisecond synaptic timing in binaural coincidence detectors by Kv1 channels"
% Paul J Mathews, Pablo E Jercog,	John Rinzel, Luisa L Scott, Nace L Golding
% Nature Neuroscience 13 601-609 (2010)

% Simulation code by Joshua H Goldwyn
% Submitted to ModelDB 1/14/13 by Joshua H Goldwyn [jgoldwyn@nyu.edu]

close all
clear all

    
%%% Set Parameters %%%
% Example: MONOLATERAL 500 Hz EXCITATION %
tEnd = 7.;             % simulation duration [ms]
stimType = 'left';     % monolateral excitation
gE = 10;               % excitatory conductance [mS / cm2]
tauE = 0.2;            % excitatory time constant (alpha function) [ms]
csynE = [2 22];        % location of excitation (compartment number)
gI = 0;                % inhibitory conductance [mS / cm2]
tauI = [.4 2];         % inhibitory time constants (double exponential function) [ms]
csynI = [12];          % location of inhibition (compartment number)
synFreq = [500 501]; % EPSP frequency (Hz) for each dendrite. inhibition freq is first entry
synDelay = [.0 .0];    % Delay of EPSP onset times in each dendrite [ms]
inhibDelay = 0;        % Delay of inhibition relative to excitation in first entry of synDelay
FreezeKLT = 0;         % Whether to Freeze KLT conductance at rest (0=No)
rB = 11;               % radius of extracellular virtual cylinder (must be larger than soma radius = 10) [micro m]

%%% Run model %%%
out = MSO_dae(tEnd, stimType, gE, tauE, csynE, gI, tauI, csynI, synFreq, synDelay, inhibDelay, FreezeKLT, rB);

%%% Surface Plot results %%%
% NOTE: Ve is not extended to ground, x-dimension is size of the neuron model
subplot(1,2,1)
surf(out.x,out.t,out.Vm), shading('flat')
set(gca,'FontSize',18)
xlabel('Distance from soma (\mum)')
ylabel('Time (ms)')
zlabel('Vm (mV)')
title('Membrane Potential','FontSize',24)

subplot(1,2,2)
surf(out.x,out.t,out.Ve), shading('flat')
set(gca,'FontSize',18)
xlabel('Distance from soma (\mum)')
ylabel('Time (ms)')
zlabel('Ve (mV)')
title('Extracellular Potential','FontSize',24)

set(gcf,'position',[25         291        1399         404])

