function trial_num = read_bit_code(x, sampleRate, bitTime, gapTime, preTime)
% 
% sampleRate in Hz.
% bitTime, gapTime, preTime all in milliseconds.
%
% DHO, 4/08
%

nbits = 10;

if length(x) > 0.08*sampleRate % Bit code happens in first 80 ms now.
    x = x(1:(0.08*sampleRate));
end

timeFromEdgeInMs = 1.5;
plusMinusTimeInMs = 2.0;


samplesFromEdge = timeFromEdgeInMs*sampleRate/1000;
plusMinusSamples = plusMinusTimeInMs*sampleRate/1000;

bitInd = zeros(1,nbits);
bitInd(1) = preTime*sampleRate/1000 + samplesFromEdge;
bitInd(2) = preTime*sampleRate/1000 + (bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(3) = preTime*sampleRate/1000 + 2*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(4) = preTime*sampleRate/1000 + 3*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(5) = preTime*sampleRate/1000 + 4*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(6) = preTime*sampleRate/1000 + 5*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(7) = preTime*sampleRate/1000 + 6*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(8) = preTime*sampleRate/1000 + 7*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(9) = preTime*sampleRate/1000 + 8*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
bitInd(10) = preTime*sampleRate/1000 + 9*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;

bitIndWindow = zeros(2*plusMinusSamples+1,length(bitInd));

for k=1:length(bitInd)
    bitIndWindow(:,k) = (bitInd(k)-plusMinusSamples):(bitInd(k)+plusMinusSamples);
end

% figure; plot(x); hold on
% plot(bitInd, repmat(5,size(bitInd)),'r*')

trial_num = LCA.binvec2dec((max(x(bitIndWindow)) > 1));











% 
% 
% bitInd = zeros(1,10);
% bitInd(1) = preTime*sampleRate/1000 + samplesFromEdge;
% bitInd(2) = preTime*sampleRate/1000 + (bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(3) = preTime*sampleRate/1000 + 2*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(4) = preTime*sampleRate/1000 + 3*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(5) = preTime*sampleRate/1000 + 4*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(6) = preTime*sampleRate/1000 + 5*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(7) = preTime*sampleRate/1000 + 6*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(8) = preTime*sampleRate/1000 + 7*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(9) = preTime*sampleRate/1000 + 8*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;
% bitInd(10) = preTime*sampleRate/1000 + 9*(bitTime+gapTime)*sampleRate/1000 + samplesFromEdge;

% trial_num = binvec2dec((x(bitInd)>1)');