% Computing spike train correlation across brain states
% As in Peyrache et al., Nat Neurosci, 2015

% Adrien Peyrache 2017 | peyrachelab.com | github.com/PeyracheLab
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

load DataHDCells_RunSleep.mat
% File contains the binned spike trains of 19 HD cells simultaneously
% recorded during RUN, REM sleep and Slow-Wave sleep (or non-REM sleep)
% Qrun/rem/sws are 2D matrices (bins X cells).
% time bins are 10ms.

%Index of HD cells coding for the same direction
%cellIx = [5 6];

%Index of HD cells coding for opposite direction
cellIx = [8 17];

[Hrun,b]    = xcorr ( Qrun(:,cellIx(1)) , Qrun(:,cellIx(2)) , 200 , 'coeff');
Hrem        = xcorr ( Qrem(:,cellIx(1)) , Qrem(:,cellIx(2)) , 200 , 'coeff');
Hsws       = xcorr ( Qsws(:,cellIx(1)) , Qsws(:,cellIx(2)) , 200 , 'coeff');

t = b*0.01; %time bins are 10ms

figure(1),clf
plot(t,Hrun,'b')
hold on, plot(t,Hrem,'r')
hold on,plot(t,Hsws,'k')
legend('RUN','REM','SWS (non-REM)','location','eastoutside')
xlabel('Time-lag (s)')
ylabel('Correlation')

