
%  WhiskerSignalTrialArray methods:
%
% [Y,X,T] = get_cached_follicle_coords(obj,tid)
% [Y,X,T] = get_follicle_coords(obj,tid,extrap_distance)
% [Y,T] = get_follicle_translation(obj,tid)
% obj = recompute_cached_follicle_coords(obj, extrap_distance, varargin)
% 
%  WhiskerSignalTrial methods:
%
% [y,x,t] = get_cached_follicle_coords(obj,tid)
% [y,x,t] = get_follicle_coords(obj,tid,extrap_distance)
% [y,t] = get_follicle_translation(obj,tid)
% obj = recompute_cached_follicle_coords(obj, extrap_distance, varargin)
% plot_follicle_position_time_projection(obj, tid, varargin)


% Possible workflow:

load('whiskerSignalTrialArray_DO79','ws') % Load a WhiskerSignalTrialArray that does not yet
                                     % have pre-computed (cached) follicle position.
% tid = 0; 
tid = [0 1 2];
extrap_distance = 5; % Extrapolate 5 pixels beyond intersection of whisker with mask.

ws = ws.recompute_cached_follicle_coords(extrap_distance, tid); % Compute and cache follicle (x,y)
                                                                % coordinates for every trial in this
                                                                % WhiskerSignalTrialArray
save('whiskerSignalTrialArray_DO79_Follicle2','ws') % Resave, now with cached follicle positions.

%-------------------------------------------------------
% Now, can get/plot translation of the follicle:
[Y,T] = ws.get_follicle_translation(tid); % Do whole WhiskerSignalTrialArray
figure; plot(T{1},Y{1},'k-'); % Plot translation for first trial in WhiskerSignalTrialArray

% Or, get the precomputed (x,y) coordinates:
[Y,X,T] = ws.get_cached_follicle_coords(tid); % or, for the whole WhiskerSignalTrialArray

%--------------------------------------------------------
load('whiskerSignalTrialArray_DO79_Follicle')% A WhiskerSignalTrialArray *with* cached follicle positions.

trajectory_nums = {[0 1 2],{'D4','D3','D2'}};  % D4=0,D3=1,D2=2
roiPos0=[1 12]; roiCurv0=[1 40]; 
roiPos1=[1 12]; roiCurv1=[1 40]; 
roiPos2=[1 12]; roiCurv2=[1 40]; 
roi = {trajectory_nums{1},[roiPos0 roiCurv0], [roiPos1 roiCurv1], [roiPos2 roiCurv2]};

x = ws.trials{10}; 

figure('Color','white'); hold on
x.plot_fitted_whisker_time_projection(2, 'b-',[],{roiPos2,'k-'})
x.plot_fitted_whisker_time_projection(1, 'r-',[],{roiPos1,'k-'})
x.plot_fitted_whisker_time_projection(0, 'g-',[],{roiPos0,'k-'})
set(gca,'XLim',[0 150],'YLim',[0 200])
x.plot_mask(0,'k-');

% Plot follicle position for whisker 0:
x.plot_follicle_position_time_projection(0, 'm.')


% Or, could "manually" get the precomputed (x,y) coordinates for one trial and plot:
[y,x,t] = x.get_cached_follicle_coords(0);  
plot(x,y,'y.')

%--------------------------------------------------------
%-- Examples of various uses of new functions:
%--------------------------------------------------------
load('whiskerSignalTrialArray_DO79') % A WhiskerSignalTrialArray without cached follicle positions.
obj = ws.trials{40};
tid = 0;
extrap_distance = 10;


[y,x,t] = obj.get_follicle_coords(tid,extrap_distance);
[Y,X,T] = ws.get_follicle_coords(tid,extrap_distance);


obj = obj.recompute_cached_follicle_coords(extrap_distance, tid); % Do one WhiskerSignalTrial
ws = ws.recompute_cached_follicle_coords(extrap_distance, tid); % Do whole WhiskerSignalTrialArray

% Get the precomputed (x,y) coordinates:
[y,x,t] = obj.get_cached_follicle_coords(tid); % for one WhiskerSignalTrial
[Y,X,T] = ws.get_cached_follicle_coords(tid); % or, for the whole WhiskerSignalTrialArray

%---
load('whiskerSignalTrialArray_DO79_Follicle')% A WhiskerSignalTrialArray *with* cached follicle positions.
obj = ws.trials{40};
tid = 0;

[y,t] = obj.get_follicle_translation(tid); % Do one WhiskerSignalTrial
[Y,T] = ws.get_follicle_translation(tid); % Do whole WhiskerSignalTrialArray

figure; plot(T{1},Y{1},'k-'); % Plot translation for first trial in WhiskerSignalTrialArray

% Get the precomputed (x,y) coordinates for one trial:
[y,x,t] = obj.get_cached_follicle_coords(tid);
