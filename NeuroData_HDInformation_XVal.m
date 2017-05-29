% Estimating HD Information by a cross-validating procedure
% Following Harris et al., Nature, 2003

% Adrien Peyrache 2017 | peyrachelab.com | github.com/PeyracheLab
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

clear all
load('DataHDCells_TuningCurves')

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
Lf = zeros(10,1);

for ep=1:10
    
    %Range of test set (offset by beginning of recording)
    tStart  = (ep-1)*dT + ang(1,1);
    tEnd    = ep*dT + ang(1,1);
    
    % Defining training set (everything except test set)
    ixAngTraining = ang(:,1)<tStart | ang(:,1)>=tEnd;
    ixSpkTraining = spk(:,1)<tStart | spk(:,1)>=tEnd;
    
    % Estimating HD tuning curve on the training set
    [hd,angBins] = hdTuningCurve(ang(ixAngTraining,2),spk(ixSpkTraining,2),6);
    
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
    
    Lf(ep) = - sum(intensityFct - firingRate)/Fs + sum(logTermF - log2(firingRate));
    
end    

Lf = sum(Lf)/T;
LfSpk = Lf/fr;

fprintf('Cross-validated information rate is %f (bit/sec)\n',Lf)
fprintf('Cross-validated information per spike is %f (bit/spk)\n\n',LfSpk)

