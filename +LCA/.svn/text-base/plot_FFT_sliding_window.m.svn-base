function [] = plot_FFT_sliding_window(f,t,YMat,fmin,fmax,FFTWindowSize)
%
%
%
%
% DHO, 3/08.
%

ind = find(f>=fmin & f<=fmax);
M = YMat(ind,:); 
f = f(ind);

F = repmat(f',1,size(M,2));
T = repmat(t,size(M,1),1);

% Shift time axis by 1/2 the bin width to center each pixel on the window
% of time in which its value was computed:
T = T+(FFTWindowSize./2);

az=0; el=90;
surf(T,F,M,'EdgeColor','none'); grid off; view(az,el);
ylim([fmin fmax])