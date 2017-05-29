function [hd,angBins] = hdTuningCurve(ang,angSpk,angBin,varargin)

% Estimating HD tuning curves
%
% USAGE
%   [hdTuning,angBins] = hdTuningCurve(ang,angSpk,angBin,Fs)
%
% Inputs:
%   ang:            a vector of head-direction from video-tracking
%   angSpk:         head-direction at times of spikes
%   angBin:         width of angular bnis (in degrees)
%   Fs (optional):  video sampling frequency (default: 39.0625)
%
% Outputs:
%   hdTuning:       tuning curve
%   angBins:        angular bins

% Adrien Peyrache 2017 | peyrachelab.com | github.com/PeyracheLab
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

%Default sampling frequency
Fs = 39.065;
if ~isempty(varargin)
    Fs = varargin{1};
end

da = angBin * pi/180; %bin width in degrees
angBins = [da/2:da:2*pi-da/2];

%Occupancy
histAng = hist(ang,angBins);

spkPerAng = hist(angSpk,angBins);
hd = spkPerAng./histAng * Fs;