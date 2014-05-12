%
%
%
%
% DHO, 5/08.
% SAH, 10/11
%
%
%
%
classdef ClustTrialArray < handle

    properties
        cellNum = '';
        shankNum = '';
        clustTrials = {};
        xsgSweepNums = [];
    end

    properties (Dependent = true)
        trialNums
        spikeRatesInHz
    end

    methods (Access = public)
        function obj = ClustTrialArray(clust_trials, cell_num, shankNum)
            %
            %   function obj = SpikesTrialArray(spikes_trials, cell_num, shankNum)
            %
            %    spikes_trials: cell array of SpikesTrial objects.
            %
            %
            if nargin > 0
                obj.cellNum = cell_num;
                obj.shankNum  = shankNum;
                obj.clustTrials = clust_trials;
            end
        end

        function r = length(obj)
            r = length(obj.clustTrials);
        end
        
        function r = peakInstRate(obj)
            if isempty(obj.clustTrials)
                 r = [];
            else
                r = cellfun(@(x) x.peakInstRate, obj.clustTrials);             
            end
        end
        
        function r = mergeAdjSpikeTimes(obj, nsamps)
            if ~isempty(obj.clustTrials)
                for k=1:length(obj)
                    obj.clustTrials{k} = obj.clustTrials{k}.mergeAdjSpikeTimes(nsamps);
                end
            end
            r = obj;
        end
    end


    methods % Dependent property methods; cannot have attributes.
        function value = get.trialNums(obj)
            if ~isempty(obj.clustTrials)
                value = cellfun(@(x) x.trialNum, obj.clustTrials);
            else
                value = [];
            end
        end
        
        function value = get.spikeRatesInHz(obj)
            if ~isempty(obj.clustTrials)
                value = cellfun(@(x) x.spikeRateInHz, obj.clustTrials);
            else
                value = [];
            end
        end

    end
end

