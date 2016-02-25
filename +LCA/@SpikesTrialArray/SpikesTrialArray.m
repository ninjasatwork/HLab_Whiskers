%
%
%
%
% DHO, 5/08.
%
%
%
%
classdef SpikesTrialArray < handle

    properties
        cellNum = '';
        cellCode = '';
        spikesTrials = {};
        xsgSweepNums = [];
    end

    properties (Dependent = true)
        trialNums
        spikeRatesInHz
    end

    methods (Access = public)
        function obj = SpikesTrialArray(spikes_trials, cell_num, cell_code)
            %
            %   function obj = SpikesTrialArray(spikes_trials, cell_num, cell_code)
            %
            %    spikes_trials: cell array of SpikesTrial objects.
            %
            %
            if nargin > 0
                obj.cellNum = cell_num;
                obj.cellCode  = cell_code;
                obj.xsgSweepNums = cellfun(@(x) x.xsgFileNum, spikes_trials);
                obj.spikesTrials = spikes_trials;
            end
        end

        function r = length(obj)
            r = length(obj.spikesTrials);
        end
        
        function r = peakInstRate(obj)
            if isempty(obj.spikesTrials)
                 r = [];
            else
                r = cellfun(@(x) x.peakInstRate, obj.spikesTrials);             
            end
        end
        
        function r = mergeAdjSpikeTimes(obj, nsamps)
            if ~isempty(obj.spikesTrials)
                for k=1:length(obj)
                    obj.spikesTrials{k} = obj.spikesTrials{k}.mergeAdjSpikeTimes(nsamps);
                end
            end
            r = obj;
        end
    end


    methods % Dependent property methods; cannot have attributes.
        function value = get.trialNums(obj)
            if ~isempty(obj.spikesTrials)
                value = cellfun(@(x) x.trialNum, obj.spikesTrials);
            else
                value = [];
            end
        end
        
        function value = get.spikeRatesInHz(obj)
            if ~isempty(obj.spikesTrials)
                value = cellfun(@(x) x.spikeRateInHz, obj.spikesTrials);
            else
                value = [];
            end
        end

    end
end

