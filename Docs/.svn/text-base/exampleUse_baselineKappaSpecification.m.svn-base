%---- Examples of different usages of measuring baseline kappa for subtraction from kappa
%     to give deltaKappa. 
d = 'C:\Documents and Settings\oconnord\My Documents\MATLAB\testOneByOneConstruction';
trajectory_nums = [0 1];

% First way:
Whisker.makeAllDirectory_WhiskerTrialLite(d,'r_in_mm',2.33,'calc_forces',true,'baseline_time_or_kappa_value',[0 .1]);
%----
% Second way:
% Or, do this another way by getting kappa from a single trial:
% Get kappa value directly for TID 0 from a single loaded WhiskerSignalTrial:
load('C:\Documents and Settings\oconnord\My Documents\MATLAB\testOneByOneConstruction\091101-JF25399b_006_WST.mat') % Loads object names 'ws'.
tid = 0; r_in_mm = 2.33; pxPerMm = 13.74; 
[thetap,kappap,y,x,t] = ws.get_theta_kappa_at_roi_point(tid,r_in_mm);
baselineKappa =  nanmean(kappap(t <= 0.1)) * pxPerMm; % **** Make sure to multiply to pxPerMm to convert 1/pixels to 1/mm,
                                                      % since kappap is in 1/pixels but 'baseline_time_or_kappa_value' must be in 1/mm.****
Whisker.makeAllDirectory_WhiskerTrialLite(d,'r_in_mm',2.33,'calc_forces',true,'baseline_time_or_kappa_value',baselineKappa);
%-----------

% View the resulting timeseries:
wl = Whisker.WhiskerTrialLiteArray(d);
tid=0; 
Whisker.view_WhiskerTrialLiteArray(wl,tid)

