classdef SiTrial < handle
    %
    %
    % Class to bind together behavior, electrophysiology, and videography
    % data.
    %
    % 
    %
    %
    % DHO, 5/08.
    %
    properties
        trialNum = [];
        behavTrial = []; % If initialized to empty Solo.BehavTrial, then save/load don't work properly. Same for clustTrial, whiskerTrial.
        shanksTrial = []; % LCA.shanksTrial;
        whiskerTrial = [];% Whisker.WhiskerSignalTrial or Whisker.WhiskerTrialLite;
    end

    properties (Dependent = true)
        cellNum
        shankNum
        mouseName
        sessionName
        trialType % 1 for go trial, 0 for no-go trial.
        trialCorrect % 1 for correct, 0 for incorrect.
        trialOutcome

        spikeRateInHz
        answerLickTime
        beamBreakTimes
        pinDescentOnsetTime
        pinAscentOnsetTime
        useFlag % 1 if useFlags for component trial objects are all one, 0 otherwise.
    end

    methods (Access = public)
        function obj = SiTrial(behav_trial, shanks_trial, varargin)
            %
            % function obj = Trial(behav_trial, clust_trial, varargin)
            %
            % behav_trial: Solo.BehavTrial object.
            % shanks_trial: LCA.ShanksTrial object.
            % varargin: Optional Whisker.WhiskerSignalTrial, Whisker.WhiskerTrialLite, 
            % or Whisker.WhiskerMeasurementsTrial object.
            %
            
            if nargin > 0

                if ~isa(behav_trial,'Solo.BehavTrial')
                    error('First argument must be a class of type Solo.BehavTrial')
                end
                if ~isa(shanks_trial,'LCA.ShanksTrial')
                    error('Second argument must be a class of type LCA.ShanksTrial')
                end
                if nargin > 2
                    if ~isa(varargin{1},'Whisker.WhiskerSignalTrial') && ...
                            ~isa(varargin{1},'Whisker.WhiskerTrialLite') && ...
                            ~isa(varargin{1},'Whisker.WhiskerMeasurementsTrial')
                        error(['Third argument must be a class of type Whisker.WhiskerSignalTrial, '...
                            'Whisker.WhiskerTrialLite, or Whisker.WhiskerMeasurementsTrial'])
                    end
                    obj.whiskerTrial = varargin{1};
                end

                obj.behavTrial = behav_trial;
                obj.shanksTrial = shanks_trial;

                if obj.behavTrial.trialNum ~= obj.shanksTrial.trialNum
                    error('behav_trial and shanks_trial arguments have different trialNum properties.')
                else
                    obj.trialNum = obj.behavTrial.trialNum;
                end
                if nargin > 2
                    if obj.whiskerTrial.trialNum ~= obj.behavTrial.trialNum
                        error('whisker_trial trialNum has trialNum mismatch.')
                    end
                end

            end
        end

        function r = spikeRateInHzTimeWindow(obj, clust, startTimeInSec, endTimeInSec)
            %
            % r = spikeRateInHzTimeWindow(cluster numbers, startTimeInSec, endTimeInSec)
            %
            % To get all clusters, pass "[]" to clust
            % Start and stop times are inclusive.
            %
            %

            if isempty(obj.shanksTrial)
                r = NaN;
                return
            end

            sampleRate = obj.shanksTrial.sampleRate;
            sweepLengthInSamples = obj.shanksTrial.sweepLengthInSamples;
             
            spikeTimesInSec = cellfun(@(x)x.spikeTimes ./ sampleRate,obj.shanksTrial.clustData,'UniformOutput',0);

            if startTimeInSec < 0 || startTimeInSec > sweepLengthInSamples / sampleRate
                error('Invalid value of parameter startTimeInSec.')
            end

            if endTimeInSec < 0 || endTimeInSec > sweepLengthInSamples / sampleRate
                error('Invalid value of parameter endTimeInSec.')
            end

            if startTimeInSec >= endTimeInSec
                error('Start time is greater or equal to end time.')
            end

            if length(spikeTimesInSec)==1
                if spikeTimesInSec==0
                    r = 0;
                    return
                end
            end

            nspikes = cellfun(@(x)length(find(x >= startTimeInSec & ...
                x <= endTimeInSec)),spikeTimesInSec,'UniformOutput',0);

            spikerate = cellfun(@(x)x / (endTimeInSec - startTimeInSec),nspikes,'UniformOutput',0);
            
            if isempty(clust)
                
                r = spikerate;
            
            elseif max(clust)>length(spikerate);
            
                error('Cluster number greater than available cluster data')
            else
                r = spikerate(clust);
            end
            

        end


        function r = spikeCountInTimeWindow(obj, clust, startTimeInSec, endTimeInSec)
            %
            % r = spikeRateInHzTimeWindow(cluster numbers, startTimeInSec, endTimeInSec)
            %
            % To get all clusters, pass "[]" to clust
            % Start and stop times are inclusive.
            %
            %

            if isempty(obj.shanksTrial)
                r = NaN;
                return
            end

            sampleRate = obj.shanksTrial.sampleRate;
            sweepLengthInSamples = obj.shanksTrial.sweepLengthInSamples;
             
            spikeTimesInSec = cellfun(@(x)x.spikeTimes ./ sampleRate,obj.shanksTrial.clustData,'UniformOutput',0);

            if startTimeInSec < 0 || startTimeInSec > sweepLengthInSamples / sampleRate
                error('Invalid value of parameter startTimeInSec.')
            end

            if endTimeInSec < 0 || endTimeInSec > sweepLengthInSamples / sampleRate
                error('Invalid value of parameter endTimeInSec.')
            end

            if startTimeInSec >= endTimeInSec
                error('Start time is greater or equal to end time.')
            end

            if length(spikeTimesInSec)==1
                if spikeTimesInSec==0
                    r = 0;
                    return
                end
            end

            nspikes = cellfun(@(x)length(find(x >= startTimeInSec & ...
                x <= endTimeInSec)),spikeTimesInSec,'UniformOutput',0);

            spikecount = cellfun(@(x)x,nspikes,'UniformOutput',0);
            
            if isempty(clust)
                
                r = spikecount;
            
            elseif max(clust)>length(spikecount);
            
                error('Cluster number greater than available cluster data')
            else
                r = spikecount(clust);
            end
            

        end
    end



    methods % Dependent property methods; cannot have attributes.

        function value = get.cellNum(obj)
            value = cellfun(@(x)x.cellNum,obj.shanksTrial.clustData);
        end

        function value = get.shankNum(obj)
            value = cellfun(@(x)x.shankNum,obj.shanksTrial.clustData);
        end

        function value = get.mouseName(obj)
            value = obj.behavTrial.mouseName;
        end

        function value = get.sessionName(obj)
            value = obj.behavTrial.sessionName;
        end

        function value = get.trialType(obj)
            value = obj.behavTrial.trialType;
        end

        function value = get.trialCorrect(obj)
            value = obj.behavTrial.trialCorrect;
        end
        
        function value = get.trialOutcome(obj)
            trialType = obj.behavTrial.trialType;
            trialCorrect = obj.behavTrial.trialCorrect;
            if trialType==1 && trialCorrect==1
                value='Hit';
            elseif trialType==1 && trialCorrect==0
                value='Miss';
            elseif trialType==0 && trialCorrect==1
                value='CorrectRejection';
            elseif trialType==0 && trialCorrect==0
                value='FalseAlarm';
            end
        end
        
        function value = get.answerLickTime(obj)
            value = obj.behavTrial.answerLickTime;
        end

        function value = get.beamBreakTimes(obj)
            value = obj.behavTrial.beamBreakTimes;
        end
        
        function value = get.pinDescentOnsetTime(obj)
            value = obj.behavTrial.pinDescentOnsetTime;
        end
        
        function value = get.pinAscentOnsetTime(obj)
            value = obj.behavTrial.pinAscentOnsetTime;
        end        

        function value = get.spikeRateInHz(obj)

            if ~isempty(obj.shanksTrial.clustData)
                value = cellfun(@(x)x.spikeRateInHz,obj.shanksTrial.clustData,'UniformOutput',0);
            else
                value = [];
            end
        end

        function value = get.useFlag(obj)
            if obj.behavTrial.useFlag==1 && obj.shanksTrial.useFlag==1 % ADD: obj.WhiskerTrial.useFlag==1
                value = 1;
            else
                value = 0;
            end
        end

    end
end

% SHOULD IMPLEMENT PLOTTING METHODS TOO.  TRIALARRAY OBJECT VIEWER WILL
% THEN JUST CALL THIS OBJECT'S PLOTTING METHOD AND WILL NOT HAVE TO ACCESS
% THIS OBJECT'S PROPERTIES.
