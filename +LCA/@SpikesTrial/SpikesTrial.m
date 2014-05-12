classdef SpikesTrial < handle
    %
    %
    %
    %
    % DHO, 5/08.
    %
    %
    %
    %
    properties
        trialNum = []; %: 108
        cellNum = ''; %'DO0037'
        cellCode = ''; % 'AAAB'
        xsgFileNum = []; %: 1
        useFlag = 1;
        sampleRate = []; %20000
        sweepLengthInSamples = [];
        spikeTimes = []; % [17x1 double]
    end

    properties (Dependent=true)
        spikeRateInHz
    end

    methods (Access = public)
        function obj = SpikesTrial(trial_num, cell_num, cell_code, xsg_file_num,...
                sample_rate, sweep_length_in_samples, spike_times)
            %
            %  function obj = SpikesTrial(trial_num, cell_num, cell_code, xsg_file_num,...
            %    sample_rate, spike_times)
            %
            %
            if nargin > 0
                obj.trialNum = trial_num;
                obj.cellNum = cell_num;
                obj.cellCode  = cell_code;
                obj.xsgFileNum = xsg_file_num;
                obj.sampleRate = sample_rate;
                obj.sweepLengthInSamples = sweep_length_in_samples;
                obj.spikeTimes = spike_times;
            end
        end

        function r = burstSpikeTimes(obj, varargin)
            % r = burstTimes(obj)
            %
            % Returns times (in samples) of locations of spikes occuring as
            %   part of bursts, defined as contigous spikes with inter-spike
            %   interval less than a threshold.
            % 
            %   Done as follows: all spikes in original spike times vector
            %   that are separated from another spike by less than the
            %   threshold interspike interval are kept, and the rest are
            %   removed.
            %
            % varargin{1} is optional instananeous spike rate in Hz to
            %   use as the criterion. Default is 100 Hz.
            % 
            
            if nargin > 1
                threshISI = 1 / varargin{1};
            else
                threshISI = 1 / 100;
            end
            
            if  isempty(obj.spikeTimes)
                r = [];
            elseif obj.spikeTimes==0
                r = 0;
            elseif length(obj.spikeTimes)==1 
                r = 0; % Can be no bursts if there's only 1 spike.
            else
                threshNumSamples = threshISI * obj.sampleRate;
                s = obj.spikeTimes';
                x = diff(s) < threshNumSamples;
                r = s(or([0 x], [x 0]))';
%                 r = s ./ obj.sampleRate;
            end
        end
        
        function r = getSpikeTimes(obj)
            %
            % r is in seconds.
            %
            if  isempty(obj.spikeTimes)
                r = [];
            elseif obj.spikeTimes==0
                r = 0;
            else
                s = obj.spikeTimes';
                r = s ./ obj.sampleRate;
            end
        end
        
        function r = mergeAdjSpikeTimes(obj, nsamps)
            %
            % Merges spike times separated by nsamps samples
            % or less, on assumption (which user must
            % be sure is correct) that these are spuriosly
            % considered separate spikes. Merges by keeping the
            % first spike time.
            %
            % E.g., if spike times are recorded on three samples: [1 0 1]
            % and you want these merged, call this function with nsamps=2.
            %
            % Returns new, corrected SpikesTrial.
            %
            if  isempty(obj.spikeTimes)
                r = obj;
                return
            elseif obj.spikeTimes==0
                r = obj;
                return
            end
            s = obj.spikeTimes;
            ds = [0; diff(s)];
            ind = ds==0 | ds>nsamps;
            obj.spikeTimes = s(ind);
            ndiscard = sum(~ind);
            if ndiscard > 0
                disp(['Found ' int2str(ndiscard) ' spikes separated by ' int2str(nsamps) '; discarded 2nd of each.'])
            end
            r = obj;
        end
        
        function r = getBurstSpikeTimes(obj)
            %
            % r is in seconds.
            %
            if  isempty(obj.spikeTimes)
                r = [];
            elseif obj.spikeTimes==0
                r = 0;
            else
                s = obj.burstSpikeTimes';
                r = s ./ obj.sampleRate;
            end
        end
        
        function r = peakInstRate(obj)
            % function r = peakInstRate()
            %
            % Peak instantaneous firing rate in Hz for the trial.
            % Equal to 1 / (minimum interspike interval) if >= 2 spikes,
            % to 1 / sweepLength (eg, 5 s) if 1 spike, or to 0 if 
            % there are 0 spikes.
            %
            if  isempty(obj.spikeTimes)
                r = [];
            elseif obj.spikeTimes==0
                r = 0;
            elseif length(obj.spikeTimes)==1 
                r = 1 / (obj.sweepLengthInSamples ./ obj.sampleRate);
            else
                st = obj.spikeTimes ./ obj.sampleRate;
                r = 1 / min(diff(st));
                if r > 500
                    beep
                end
            end
        end
              
        function plot_trial_events(obj)
            cla
            ymin = .5; ymax = 1.5;
            st = obj.spikeTimes ./ obj.sampleRate;
            sweepEndTime = obj.sweepLengthInSamples ./ obj.sampleRate;
            plot(st, 1, 'ko')
            xlim([0 sweepEndTime])
            ylim([ymin ymax])
            set(gca,'YTick',[],'box','off','TickDir','out')
            title(['TrialNum=' int2str(obj.trialNum)])
        end
       
    end

    methods % Dependent property methods; cannot have attributes.

        function value = get.spikeRateInHz(obj)
            if  isempty(obj.spikeTimes)
                value = [];
            elseif obj.spikeTimes==0
                value = 0;
            else
                sweepLengthInSec = obj.sweepLengthInSamples ./ obj.sampleRate;
                value = length(obj.spikeTimes) ./ sweepLengthInSec;
            end
        end

    end
end

