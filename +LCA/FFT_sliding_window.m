function [f,t,YMat] = FFT_sliding_window(y, Fs, windowSize, stepSize)
%
% y:  the signal to compute FFT on.
% Fs: sampling rate in Hz.
% windowSize: time of window in which to compute FFT in seconds.
% stepSize: how far to move in seconds before each successive FFT.
%
%
% DHO, 3/08
%

windowSize = windowSize*Fs;
stepSize = stepSize*Fs;

totalLength = length(y);
nSteps = (totalLength-windowSize)/stepSize;

NFFT = 2^nextpow2(windowSize);
h=hamming(windowSize,'periodic');


YMat = zeros(NFFT/2,nSteps);

% Mean-subtract signal to eliminate DC term of transform:
n=1;
for k=1:nSteps
    yy = y(n:(n+windowSize-1));
    yy = yy-mean(yy);

    % % Multiply by Hamming window to reduce spectral leakage:
    yy = yy.*h;

    Y = fft(yy,NFFT);
    Y = 2*abs(Y(1:NFFT/2));

    YMat(:,k) = Y';
    
    n = n + stepSize;
end

f = Fs/2*linspace(0,1,NFFT/2);  % frequency in Hz.
t = stepSize*(0:(nSteps-1))/Fs; % time in seconds.

f = fliplr(f);
YMat = flipud(YMat);