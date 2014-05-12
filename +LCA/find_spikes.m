function r = find_spikes(x,spikeThreshold, artifactThreshold)
% 
%
% DHO, 12/07
%


r = [0]; % returns this if finds no spikes

% [lmval,indd] = lmax(x); 

maximaInds = findmaxima(x);
maximaVals = x(maximaInds);

tmp = find(maximaVals>spikeThreshold);


spikeInds = maximaInds(tmp);

minimaInds = findminima(x);
minimaInds = minimaInds(x(minimaInds)<artifactThreshold);
artifacts = {};
for i = 1:length(minimaInds)
    artifacts{i} = spikeInds(abs(spikeInds-minimaInds(i))<25)';
end

spikeInds = setdiff(spikeInds,[artifacts{:}]);

% spikeInds = indd(find(lmval>thresh));

if ~isempty(spikeInds)
    r = spikeInds;
end




function maxima = findmaxima(x)
%FINDMAXIMA  Find location of local maxima
%  From David Sampson
%  See also FINDMINIMA

% Unwrap to vector
x = x(:);
% Identify whether signal is rising or falling
%upordown = sign(smooth(diff(x),3));
upordown = sign(diff(x));
% Find points where signal is rising before, falling after
maxflags = [upordown(1)<0; diff(upordown)<0; upordown(end)>0];
maxima   = find(maxflags);


function minima = findminima(x)
%FINDMAXIMA  Find location of local maxima
%  From David Sampson
%  See also FINDMINIMA

% Unwrap to vector
x = x(:);
% Identify whether signal is rising or falling
upordown = sign(smooth(diff(x),3));
% Find points where signal is rising before, falling after
minflags = [upordown(1)>0; diff(upordown)>0; upordown(end)<0];
minima   = find(minflags);