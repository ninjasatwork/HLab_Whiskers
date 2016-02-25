function [t, sampleRate, varargout] = load_xsg_trace(cellnum,code,tracenum,chans)
%
%
%
% t: time in minutes from 1-Jan-0000 reference (datenum(datevec(timeOfTrigger))*60*24).
% varargout: voltage values from the channels specified in 'chans'.
%
% 12/07, DHO
%

if tracenum < 10
    filler = '000';
elseif tracenum <100
    filler = '00';
elseif tracenum <1000
    filler = '0';
end

d = load([cellnum filesep cellnum code filler int2str(tracenum) '.xsg'],'-mat');

t = datenum(d.header.ephys.ephys.triggerTime)*60*24;
% t = datenum(d.header.ephys.ephys.triggerTime)*60*60*24;

sampleRate = d.header.ephys.ephys.sampleRate;

if length(chans)==1
    if chans==1 
        varargout{1} = d.data.ephys.trace_1;
    elseif chans==2
        varargout{1} = d.data.acquirer.trace_1;
    elseif chans==3
        varargout{1} = d.data.ephys.trace_3;
    end
else
    if max(chans)>3
        error('Invalid channel number')
    end

       eval(['varargout{1} = d.data.ephys.trace_1 ;'])
       eval(['varargout{2} = d.data.acquirer.trace_3 ;'])
end
        

