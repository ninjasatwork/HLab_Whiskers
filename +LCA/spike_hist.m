function [n,x] = spike_hist(spikeTimes, nbins)
%
%
%
%
st = [];

for k=1:length(spikeTimes)
    if ~isempty(spikeTimes{k})
        st = [st; spikeTimes{k}];
    end
end

[n,x] = hist(st, nbins);
