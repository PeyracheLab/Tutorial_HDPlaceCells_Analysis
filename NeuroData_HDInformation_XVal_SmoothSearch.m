% Cross-Validated Information measure
% 'Grid-search' of optimal angle smoothing window

% Adrien Peyrache 2017
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

clear all

load('DataHDCells_TuningCurves')

binSize = 1; %in degrees
%Smoothing windows for the tuning curve (in degrees)
smoothWin = [1 1 3 6 12 24 48];

Fs = 39.0625; % video-tracking sampling frequency
% length of recording
T = size(ang,1)/Fs;

%Let's look at the information of a cell
cellIx = 5;
spk = spikeTA{cellIx};

% Number of spikes
N = size(spk,1);
% Average firing rate
fr = N/T;

% We want to do a 10-fold cross-validation
% Length of 1/10th of the epoch
dT = T/10;

%Likelihood function is divided in 10 values
Lf = zeros(10,length(smoothWin));
for lambda = 1:length(smoothWin)
for ep=1:10
    
    %Range of test set (offset by beginning of recording)
    tStart  = (ep-1)*dT + ang(1,1);
    tEnd    = ep*dT + ang(1,1);
    
    % Defining training set (everything except test set)
    ixAngTraining = ang(:,1)<tStart | ang(:,1)>=tEnd;
    ixSpkTraining = spk(:,1)<tStart | spk(:,1)>=tEnd;
    
    % Estimating HD tuning curve on the training set
    [hd,angBins] = hdTuningCurve(ang(ixAngTraining,2),spk(ixSpkTraining,2),binSize);
    
    %Smoothing the tuning window
    l       = length(hd);
    hdTmp   = [fliplr(hd) hd fliplr(hd)];
    Npoints = round(10*smoothWin(lambda)/binSize);
    gw      = gausswin(Npoints,5); %alpha = 0.05
    gw      = gw/sum(gw);
    hdTmp   = convn(hdTmp(:),gw,'same');
    hd      = hdTmp(l+1:2*l);
    
    % Defining test set
    ixAngTest = ang(:,1)>=tStart & ang(:,1)<tEnd;    
    ixSpkTest = spk(:,1)>=tStart & spk(:,1)<tEnd;
        
    % Computing the expected firing on the test set
    angTest = ang(ixAngTest,2);
    
    % index of angBins corresponding to each angle during the test
    normAng = (angTest - min(angTest))/(max(angTest) - min(angTest));
    xx = floor((length(angBins)-1)*normAng)+1;
    
    % expected firing rate (intensity function)
    intensityFct = hd(xx);
    
    % Evaluation of the intensity function at spike times
    % First, index of 'ang' corresponding to spikes
    tAngTest = ang(ixAngTest,1);
    tSpkTest = spk(ixSpkTest,1);
    spkAngIx = zeros(length(tSpkTest),1);

    for tt=1:size(tSpkTest,1)
        [~,spkAngIx(tt)] = min(abs(tAngTest(:,1)-tSpkTest(tt)));
    end
    
    % Likelihood function
    % 'Flat' firing rate
    firingRate = length(tSpkTest)/dT;
    
    % Log terms for the intensity functions
    logTermF = log2(intensityFct(spkAngIx));
    logTermF(intensityFct(spkAngIx)==0) = 0;
    
    Lf(ep,lambda) = - sum(intensityFct - firingRate)/Fs + sum(logTermF - log2(firingRate));
    
end    
end

Lf = nansum(Lf)/T;
LfSpk = Lf/fr;

figure(1),clf
semilogx(smoothWin,LfSpk)
xlabel('Smoothing Window (degrees)')
ylabel('Cross-Validated Information (bit per spk)')
title('Cross-Validated Information')