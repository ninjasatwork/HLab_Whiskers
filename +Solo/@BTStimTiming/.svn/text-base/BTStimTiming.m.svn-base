%
%
%
%
% DHO, 8/10.
%
%
%
%
classdef BTStimTiming < handle

    properties

        mouseName = '';
        sessionName = '';
        trialNum  = [];
        trialType = ''; %  {'VirtUnpunishedNogo', 'VirtUnrewardedGo','Unpunished Nogo','Unrewarded Go','Nogo','Go'}
        response = []; % lick: 1, nolick: 0.
        
        trialEvents = [];
        trialEmbCLogItems = {};
        
        nextTrialEvents = []; % Need to store trialEvents for next trial because licks in state 35 will largely
                              % be associated with next trial.
        useFlag   = 1; % set to 0 to mark bad trials.
        sessionType = '';
        extraITIOnErrorSetting = [];
        samplingPeriodTimeSetting = [];
        waterValveTimeSetting = [];
        answerPeriodTimeSetting = []; % 2 sec minus samplingPeriodTimeSetting
        motorPosition = [];
        nogoPosition = [];
        goPosition = [];      
    end


    properties (Dependent = true, SetAccess = private)
        beamBreakTimes = [];    % Note that the last trial acquired in Solo will lack beambreaks
                                % occuring during the subsequent inter-trial-interval (i.e. state 35).
        trialTriggerTimeEPHUS = [];
        trialTriggerTimeCamera = []; 
        rewardTime = []; % [startTime stopTime]
        airpuffTimes = {}; % Cell array of {[startTime1 stopTime1],[startTime2 stopTime2],...} 
        pinMoveOutOfReachOnsetTime = []; ;
        pinMoveIntoReachOnsetTime = []; 
        
        samplingPeriodTime = []; % Retrieved from event matrix
        answerPeriodTime = []; % Retrieved from event matrix; maximum of 
                                % answerPeriodTimeSetting, but ends when mouse licks.
                                % Gives reaction time on Go trials from end
                                % of sampling/grace period.

        answerLickTime = []; % Empty if a miss or correct rejection.
        trialStartTime = [];
        drinkingTime = []; % 2 s minus water valve time, to give mouse time to drink before proceeding w/ next trial.
        timeoutPeriodTimes = [];
         
        targetCrossingTimes = []; % from EmbC  varLog
        stimTimes473nm = []; % from EmbC varLog
        
        trialCorrect 
        stimGiven
        numStims
%         stimGivenButOnlyAfterAnswer
        
    end

    methods (Access = public)

        
        function obj = BTStimTiming(mouse_name, session_name, trial_num, trial_type,...
                response, trial_events, next_trial_events, embc_log_items, varargin)
            %         function obj = BTStimTiming(mouse_name, session_name, trial_num, trial_type,...
            %                 trial_response, trial_events, next_trial_events, varargin)
            %
            %              Optional parameter/value pairs: useFlag, sessionType, extraITIOnErrorSetting,
            %              samplingPeriodTimeSetting, waterValveTimeSetting,
            %              motorPosition, nogoPosition, goPosition
            %
            if nargin == 0
                return
            end
            
            p = inputParser;
            
            p.addRequired('mouse_name', @ischar);
            p.addRequired('session_name', @ischar);
            p.addRequired('trial_num', @isnumeric);
            p.addRequired('trial_type', @ischar);
            p.addRequired('response', @isnumeric);
            p.addRequired('trial_events', @isnumeric);
            p.addRequired('next_trial_events', @isnumeric);
            p.addRequired('embc_log_items', @iscell);
            
            p.addParamValue('useFlag', 1, @isnumeric);
            p.addParamValue('sessionType', '', @ischar);
            p.addParamValue('extraITIOnErrorSetting', NaN, @isnumeric);
            p.addParamValue('samplingPeriodTimeSetting', NaN, @isnumeric);
            
            p.addParamValue('waterValveTimeSetting', NaN, @isnumeric);
            p.addParamValue('motorPosition', NaN, @isnumeric);
            p.addParamValue('goPosition', NaN, @isnumeric);
            p.addParamValue('nogoPosition', NaN, @isnumeric);
            
            p.parse(mouse_name, session_name, trial_num, trial_type,...
                response, trial_events, next_trial_events, embc_log_items,varargin{:});
            
%             disp 'List of all arguments:'
%             disp(p.Results)
            
            obj.mouseName = p.Results.mouse_name;
            obj.sessionName = p.Results.session_name;
            obj.trialNum  = p.Results.trial_num;
            obj.trialType  = p.Results.trial_type;
            obj.response = p.Results.response; % 1 for lick, 0 for no-lick.
            obj.trialEvents = p.Results.trial_events;
            obj.nextTrialEvents = p.Results.next_trial_events;   
            obj.trialEmbCLogItems = p.Results.embc_log_items;
            
            obj.useFlag   = p.Results.useFlag; % set to 0 to mark bad trials.
            obj.sessionType = p.Results.sessionType;
            obj.extraITIOnErrorSetting = p.Results.extraITIOnErrorSetting;
            obj.samplingPeriodTimeSetting = p.Results.samplingPeriodTimeSetting;
            
            obj.answerPeriodTimeSetting = 2 - obj.samplingPeriodTimeSetting;
            
            obj.waterValveTimeSetting = p.Results.waterValveTimeSetting;
            obj.motorPosition = p.Results.motorPosition;
            obj.nogoPosition = p.Results.nogoPosition;
            obj.goPosition = p.Results.goPosition;
        end
                   
        function plot_trial_events(obj)
            cla
            ymin = 0; ymax = 7; lw = 5; barw=.25;
            
            x = [obj.pinMoveOutOfReachOnsetTime, obj.pinMoveIntoReachOnsetTime]; 
            if numel(x) < 2
                return % If Solo was stopped/started sometimes the pinMoveIntoReachOnsetTime is missing
            end
%             y = 6*ones(size(x));
%             plot(x, y, 'k-','LineWidth',lw); hold on            
            xx=[x(1) x(2) x(2) x(1)]; yy= 5 + [-barw -barw barw barw]; hold on
            patch(xx,yy, 0.6*[1 1 1],'FaceAlpha',0.5,'LineStyle','none');

%             x = obj.samplingPeriodTime; 
% %             y = 5*ones(size(x));
% %             plot(x, y, 'c-','LineWidth',lw); hold on
%             xx=[x(1) x(2) x(2) x(1)]; yy= 5 + [-barw -barw barw barw];
%             patch(xx,yy, 'c','LineStyle','none');
            
            x = obj.answerPeriodTime; 
%             y = 4*ones(size(x));
%             plot(x, y, 'g-','LineWidth',lw); hold on
            xx=[x(1) x(2) x(2) x(1)]; yy= 4 + [-barw -barw barw barw]; 
            patch(xx,yy, 'k','LineStyle','none');
            
            x = obj.beamBreakTimes; 
            y = 3*ones(size(x));
            plot(x, y, 'mo','MarkerSize',9)
            
            if ~isempty(obj.rewardTime)
                x = [obj.rewardTime(1), obj.rewardTime(2)]; 
%                 y = 2*ones(size(x));
%                 plot(x, y, 'b-','LineWidth',lw)
                xx=[x(1) x(2) x(2) x(1)]; yy= 2 + [-barw -barw barw barw]; 
                patch(xx,yy, [0 96 255]./255,'FaceAlpha',1,'LineStyle','none');
            end
            
            p = obj.airpuffTimes;
            if ~isempty(p)
                for k=1:length(p)
                    puff_times = p{k};
                    x = puff_times; 
%                     y = ones(size(x));
%                     plot(x, y, 'r-','LineWidth',lw)
                    xx=[x(1) x(2) x(2) x(1)]; yy= 1 + [-barw -barw barw barw]; 
                    patch(xx,yy, [255 159 0]./255,'FaceAlpha',1,'LineStyle','none');
                end
            end
            ylim([ymin ymax]);
            
            if (obj.pinMoveIntoReachOnsetTime) > 5
                xlim([0, obj.pinMoveIntoReachOnsetTime + 0.1])
            else
                xlim([0 5])
            end
            
%             xlm = get(gca, 'XLim');
%             if xlm(2) < 5
%                 xlim([0 5])
%             end
            
            set(gca, 'YTick', 1:5,'YTickLabel', {'Airpuff','Water valve','Beam breaks', 'Answer period','Pin valve'},...
                'FontSize', 12, 'TickDir','out','Box','off')
            xlabel('Sec','FontSize',12)
            set(gcf,'Color','white')
            
            
            if obj.trialType==1
                trial_type_string = 'Go';
            else
                trial_type_string = 'Nogo';
            end
            
            if obj.trialResponse==1
                score_string = 'Licked';
            else
                score_string = 'No lick';
            end

            title(['TrialNum=' int2str(obj.trialNum) ...
                ', ' trial_type_string ', ' score_string])
               
        end
    end
    
    methods % Dependent property methods; cannot have attributes.    
        
        function value = get.trialCorrect(obj)
            tt = obj.trialType;
            lick = obj.response;
            
            go_types = {'Go','GoStim_0ms','GoStim_5ms','GoStim_10ms','GoStim_20ms','GoStim_50ms'};
            nogo_types = {'Nogo','NogoStim_0ms','NogoStim_5ms','NogoStim_10ms','NogoStim_20ms','NogoStim_50ms'};
            
            if ismember(tt,go_types) && lick==1
                value = 1;
            elseif ismember(tt,go_types) && lick==0
                value = 0;
            elseif ismember(tt,nogo_types) && lick==1
                value = 0;
            elseif ismember(tt,nogo_types) && lick==0
                value = 1;
            else
                error('Unrecognized trial type.')
            end
        end
        
%         function value = get.stimGiven(obj)
%             embc = obj.trialEmbCLogItems;
%             if sum(ismember(embc(:,2),'stim_473nm')) == 1
%                 value = true;
%             else
%                 value = false;
%             end
%         end
             
        function value = get.stimGiven(obj)
            embc = obj.trialEmbCLogItems;
            if ismember('stim_473nm',embc(:,2))
                value = true;
            else
                value = false;
            end
        end

        function value = get.numStims(obj)
            embc = obj.trialEmbCLogItems;
            value = numel(find(ismember(embc(:,2),'stim_473nm')));
        end
        
%         
%         function value = get.stimGivenButOnlyAfterAnswer(obj) % This will be obsolete for sessions after 27aug10, when stims only give before answer lick.
%             embc = obj.trialEmbCLogItems;
%             reward = obj.rewardTime;
%             puff = obj.airpuffTimes;
%             
%             if isempty(reward) & isempty(reward
%             
%             if ismember('stim_473nm',embc(:,2))
%                 
%                 value = true;
%             else
%                 value = false;
%             end
%         end
        
        function value = get.beamBreakTimes(obj)
            trialEntryInd = find(obj.trialEvents(:,1)==40,1,'first');
            breakInd = find(obj.trialEvents(:,2)==1); % Limit to events occurring after the entry to
            % to the current trial (state 40).
            breakInd = breakInd(breakInd >= trialEntryInd);
            breakTimes = obj.trialEvents(breakInd, 3) - obj.trialStartTime;
            
            % Add beam breaks occuring in the subsequent intertrial interval:
            if ~isempty(obj.nextTrialEvents)
                ITIOverInd = find(obj.nextTrialEvents(:,1)==40,1,'first');
                breakIndITI = find(obj.nextTrialEvents(:,2)==1); % Get all beam breaks assigned to next trial
                % And then restrict to those occuring before state 40 of next trial.
                % These are in the intertrial interval and we'll include them for the
                % present trial. Thus, there will be no beam breaks with negative
                % times when aligned on state 40 entry.
                breakIndITI = breakIndITI(breakIndITI < ITIOverInd); 
                breakTimesITI = obj.nextTrialEvents(breakIndITI, 3) - obj.trialStartTime;
                value = [breakTimes; breakTimesITI];
            else
                value = breakTimes;
            end
        end
        
        function value = get.rewardTime(obj)% State 43 entries and exits
            trial_events = obj.trialEvents;
            
            rowIndStart = find(trial_events(:,1)==43 & trial_events(:,2)==0, 1, 'first');
            rowIndStop = find(trial_events(:,1)==43 & trial_events(:,2)==3, 1, 'first'); % Timeout code = 3;

            if ~isempty(rowIndStart)
                value = [trial_events(rowIndStart, 3), trial_events(rowIndStop, 3)] - obj.trialStartTime;
            else
                value = [];
            end
        end
        
        function value = get.drinkingTime(obj)% State 44 entries and exits
            trial_events = obj.trialEvents;
            
            rowIndStart = find(trial_events(:,1)==44 & trial_events(:,2)==0, 1, 'first');
            rowIndStop = find(trial_events(:,1)==44 & trial_events(:,2)==3, 1, 'first'); % Timeout code = 3;

            if ~isempty(rowIndStart)
                value = [trial_events(rowIndStart, 3), trial_events(rowIndStop, 3)] - obj.trialStartTime;
            else
                value = [];
            end
        end
        
        function value = get.airpuffTimes(obj) % State 49 entries and exits.
            trial_events = obj.trialEvents;
            
            % Find first airpuff state entry, then first exit and pair them.  
            % Eliminate airpuff state event matrix entries (eg, licks in and out) 
            % already paired an any entries between paired entry/exit. 
            % Repeat until none are left.
            airpuff_events = trial_events(trial_events(:,1)==49,:);
            if isempty(airpuff_events)
                value = {};
            else
                num_exits = length(find(airpuff_events(:,2)==3));
                value = cell(1,num_exits);
                for k=1:num_exits
                    entry_ind = find(airpuff_events(:,2)==0, 1, 'first');
                    exit_ind = find(airpuff_events(:,2)==3, 1, 'first');
                    
                    entry_time = airpuff_events(entry_ind, 3);
                    exit_time = airpuff_events(exit_ind, 3);

                    value{k} = [entry_time, exit_time] - obj.trialStartTime;
                    airpuff_events = airpuff_events((exit_ind+1):end, :);
                end
            end
        end
        
        function value = get.timeoutPeriodTimes(obj) % State 45 entries and exits.
            trial_events = obj.trialEvents;
            
            % Find first timeout period state entry, then first exit and pair them.  
            % Eliminate timeout period state event matrix entries (eg, licks in and out) 
            % already paired an any entries between paired entry/exit. 
            % Repeat until none are left.
            timeout_period_events = trial_events(trial_events(:,1)==45,:);
            if isempty(timeout_period_events)
                value = {};
            else
                num_exits = length(find(timeout_period_events(:,2)==3));
                value = cell(1,num_exits);
                for k=1:num_exits
                    entry_ind = find(timeout_period_events(:,2)==0, 1, 'first');
                    exit_ind = find(timeout_period_events(:,2)==3, 1, 'first');
                    
                    entry_time = timeout_period_events(entry_ind, 3);
                    exit_time = timeout_period_events(exit_ind, 3);

                    value{k} = [entry_time, exit_time] - obj.trialStartTime;
                    timeout_period_events = timeout_period_events((exit_ind+1):end, :);
                end
            end
        end

        function value = get.samplingPeriodTime(obj)% State 41 entries and exits
            trial_events = obj.trialEvents;
            rowIndStart = find(trial_events(:,1)==41 & trial_events(:,2)==0, 1, 'first');
            rowIndStop = find(trial_events(:,1)==41 & trial_events(:,2)==3, 1, 'first'); % Timeout code = 3;
            if ~isempty(rowIndStart)
                value = [trial_events(rowIndStart, 3), trial_events(rowIndStop, 3)] - obj.trialStartTime;
            else
                value = [];
            end
        end
        
        function value = get.answerPeriodTime(obj)% State 42 entries and exits
            trial_events = obj.trialEvents;
            rowIndStart = find(trial_events(:,1)==42 & trial_events(:,2)==0, 1, 'first');
            rowIndStop = find(trial_events(:,1)==42 & ismember(trial_events(:,2), [1 2 3]), 1, 'first'); % Can exit via timeout, lick in, or lick out
            if ~isempty(rowIndStart)
                value = [trial_events(rowIndStart, 3), trial_events(rowIndStop, 3)] - obj.trialStartTime;
            else
                value = [];
            end
        end


        function value = get.pinMoveOutOfReachOnsetTime(obj) % State 41 entry
            rowInd = find(obj.trialEvents(:,1)==41,1);
            if ~isempty(rowInd)
                value = obj.trialEvents(rowInd, 3) - obj.trialStartTime;
            else
                value = [];
            end
        end
          
        function value = get.pinMoveIntoReachOnsetTime(obj) % State 48 entry
            rowInd = find(obj.trialEvents(:,1)==48,1);
            if ~isempty(rowInd)
                value = obj.trialEvents(rowInd, 3) - obj.trialStartTime;
            else
                value = [];
            end
        end

        function value = get.trialTriggerTimeEPHUS(obj) % State 40 entry
            value = 0; %obj.trialStartTime;
%             rowInd = find(obj.trialEvents(:,1)==40,1);
%             if ~isempty(rowInd)
%                 value = obj.trialEvents(rowInd, 3);
%             else
%                 value = [];
%             end
        end

        function value = get.trialTriggerTimeCamera(obj) % State 40 entry
             value = 0; %obj.trialStartTime;
%             rowInd = find(obj.trialEvents(:,1)==40,1);
%             if ~isempty(rowInd)
%                 value = obj.trialEvents(rowInd, 3);
%             else
%                 value = [];
%             end
        end
        
        function value = get.answerLickTime(obj) 
            if obj.trialType==1 && obj.trialCorrect==1 % Hit 
                value = obj.rewardTime;  % differs from reward onset time by at most 1/6000 sec (period of RTLinux server).  
                if length(value) > 1 % Can be empty in rare trials (due to stopping/starting Solo) that will be excluded later when making BehavTrialArray.
                    value = value(1);
                end
%                 value = obj.rewardTime(1); % differs from reward onset time by at most 1/6000 sec (period of RTLinux server).
            elseif obj.trialType==0 && obj.trialCorrect==0 % False Alarm.
                value = obj.airpuffTimes{1}(1); % differs from onset of first airpuff time by at most 1/6000 sec (period of RTLinux server).
            else
                value = []; % leave empty if trial is a correct rejection or a miss.
            end        
        end

        function value = get.trialStartTime(obj) % State 40 entry
            rowInd = find(obj.trialEvents(:,1)==40,1);
            if ~isempty(rowInd)
                value = obj.trialEvents(rowInd, 3);
            else
                value = [];
            end
        end

    end
end

