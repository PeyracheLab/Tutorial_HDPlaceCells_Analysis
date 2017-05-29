% Example script showing how to compute the head-direction tuning curve 

% Adrien Peyrache 2017 | peyrachelab.com | github.com/PeyracheLab
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

clear all
load('DataHDCells_TuningCurves')
% This file contains:
% spikeTA: a cell array that contains 19 elements (19 HD cells
% simultaneously recorded)
% 
% each element of the array is a 2 column matrix: spike times (in seconds)
% and head directions (in radians) at time of the spikes.
%
% ang: head-direction, a 2 column matrix of [times angle (rad)]
% X: x position, a 2 column matrix [times x-pos (cms)]
% Y: y position, a 2 column matrix [times y-pos (cms)]

Fs = 39.0625; % video-tracking sampling frequency

%angular bins
da = pi/30; %6 degrees
angBins = [da/2:da:2*pi-da/2];

%Occupancy
histAng = hist(ang(:,2),angBins);

%Let's look at the tuning curve of one cell
cellIx = 5;
spk = spikeTA{cellIx};

% histogram of the number of times the cell fired in each bin of
% head-direction
spkPerAng = hist(spk(:,2),angBins);

% now compute the tuning curve:
hdTuning = spkPerAng./histAng * Fs;

figure(1),clf
set(gcf,'Position',[62   319   783   281])
subplot(1,3,1)
    polar(angBins,spkPerAng)
    title('Number of spikes')
subplot(1,3,2)
    polar(angBins,histAng)
    title('Occupancy')
subplot(1,3,3)
    polar(angBins,hdTuning)
    title('Tuning Curve (Hz)')

