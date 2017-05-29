% Estimating HD Information from tuning curves
% Following Skaggs et al., 1993

% Adrien Peyrache 2017 | peyrachelab.com | github.com/PeyracheLab
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or

% (at your option) any later version.
clear all
load('DataHDCells_TuningCurves');
% See NeuroData_HDTuningCurve.m for mor information on the content of this file

Fs = 39.0625; % video-tracking sampling frequency
% length of recording
T = size(ang,1)/Fs;

%angular bins
da = pi/30; %6 degrees
angBins = [da/2:da:2*pi-da/2];

%Occupancy
histAng = hist(ang(:,2),angBins);

%Let's look at the information of a cell
cellIx = 5;
spk = spikeTA{cellIx};

spkPerAng = hist(spk(:,2),angBins);
hdTuning = spkPerAng./histAng * Fs;

% Number of spikes
N = size(spk,1);
% Average firing rate
fr = N/T;

% probability of occupancy:
Px = histAng./sum(histAng);

logTerm = log2(hdTuning/fr);
% Correct for undefined values
logTerm(hdTuning==0) = 0;

% Little trick to express a sum as a dot product
I = hdTuning * (logTerm.*Px)' ;

% Divide by firing rate to obtain information per spike 
Ispk = I/fr;

fprintf('Information rate is %f (bit/sec)\n',I)
fprintf('Information per spike is %f (bit/spk)\n\n',Ispk)

