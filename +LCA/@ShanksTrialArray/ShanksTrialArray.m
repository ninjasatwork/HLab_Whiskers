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
classdef ShanksTrialArray < handle

    properties
        animal = '';
        session = '';
        shankTrials = {};
    end

    properties (Dependent = true)
        trialNums
        shankNums
    end

    methods (Access = public)
        function obj = ShanksTrialArray(animal, session, shankTrials)
            %
            %   function obj = SpikesTrialArray(spikes_trials, cell_num, shankNum)
            %
            %    spikes_trials: cell array of SpikesTrial objects.
            %
            %
            if nargin > 0
                obj.animal = animal;
                obj.session  = session;
                obj.shankTrials = shankTrials;
            end
        end

        function r = length(obj)
            r = length(obj.shankTrials);
        end
        
        function r = peakInstRate(obj)
            if isempty(obj.shankTrials)
                 r = [];
            else
                r = cellfun(@(x) x.peakInstRate, obj.shankTrials);             
            end
        end
        
        function r = mergeAdjSpikeTimes(obj, nsamps)
            if ~isempty(obj.shankTrials)
                for k=1:length(obj)
                    obj.shankTrials{k} = obj.shankTrials{k}.mergeAdjSpikeTimes(nsamps);
                end
            end
            r = obj;
        end
    end


    methods % Dependent property methods; cannot have attributes.
        function value = get.trialNums(obj)
            if ~isempty(obj.shankTrials)
                value = cellfun(@(x) x.trialNum, obj.shankTrials);
            else
                value = [];
            end
        end
      
         function value = get.shankNums(obj)
            if ~isempty(obj.shankTrials)
                value = obj.shankTrials{1}.shanks;
            else
                value = [];
            end
        end

    end
end

