
%% Individual WhiskerTrial and WhiskerSignalTrial example:

clear; trajectory_nums = 0; 
fn = {'JF4793_090307_S0_trial125-cropContra'...  
'JF4793_090307_S0_trial158-cropContra'...  
'JF4793_090307_S0_trial175-cropContra'...  
'JF4793_090307_S0_trial189-cropContra'...  
'JF4793_090307_S0_trial205-cropContra'...  
'JF4793_090307_S0_trial22-cropContra'...   
'JF4793_090307_S0_trial49-cropContra'...   
'JF4793_090307_S0_trial63-cropContra',trajectory_nums};

w = Whisker.WhiskerTrialArray(fn); w.mouseName = 'JF4793'; w.sessionName = 'JF4793_090307';
% save('C:\dan\analysis\behav\analysis\behaviorPerSe\psychCurvesData\whiskerTrialArrays\whiskerTrialArray_JF4793_090307.mat','w','-v7') 
% load('C:\dan\analysis\behav\analysis\behaviorPerSe\psychCurvesData\whiskerTrialArrays\whiskerTrialArray_JF4793_090307.mat','w') 
xroiCurv=[118 300]; xroiPos=[118 159];
ws = Whisker.WhiskerSignalTrialArray(w,xroiCurv,xroiPos);
% save('C:\dan\analysis\behav\analysis\behaviorPerSe\psychCurvesData\whiskerSignalTrialArrays\whiskerSignalTrialArray_JF4793_090307.mat','ws') 
% load('C:\dan\analysis\behav\analysis\behaviorPerSe\psychCurvesData\whiskerSignalTrialArrays\whiskerSignalTrialArray_JF4793_090307.mat','ws') 

%-------------------------------------------------------------------------
p='C:\dan\analysis\behav\analysis\behaviorPerSe\analysisScripts\TrackedExampleTrials\';
% Trial is a correct rejection. 500 FPS.
fn = 'C:\dan\analysis\behav\analysis\behaviorPerSe\analysisScripts\TrackedExampleTrials\JF4004_091707_S0_trial90';


trialNum=90; mouseName='JF4004'; sessionName='091707';
trajectory_nums = [0 1 2 -1 -2 -3 ];
% w = Whisker.WhiskerTrial(fn, trialNum, trajectory_nums, mouseName, sessionName); 
% save('C:\dan\analysis\behav\analysis\behaviorPerSe\whiskerTrials\whiskerTrial_JF4004_091707_S0_trial90.mat','w','-v7')
load('C:\dan\analysis\behav\analysis\behaviorPerSe\whiskerTrials\whiskerTrial_JF4004_091707_S0_trial90.mat','w')

xroiPos0=[260 300]; xroiCurv0=[200 300]; % not using curvature from this trial.
xroiPos1=[280 320]; xroiCurv1=[250 320]; 
xroiPos2=[290 325]; xroiCurv2=[250 325]; 
xroiPosMinus1=[538 580]; xroiCurvMinus1=[538 660]; % not using curvature from this trial.
xroiPosMinus2=[538 580]; xroiCurvMinus2=[538 660]; 
xroiPosMinus3=[538 580]; xroiCurvMinus3=[538 660]; 

xroi = {trajectory_nums,[xroiPos0 xroiCurv0], [xroiPos1 xroiCurv1], [xroiPos2 xroiCurv2],...
                           [xroiPosMinus1 xroiCurvMinus1], [xroiPosMinus2 xroiCurvMinus2], [xroiPosMinus3 xroiCurvMinus3]};

ws = Whisker.WhiskerSignalTrial(w,xroi);
% save('C:\dan\analysis\behav\analysis\behaviorPerSe\whiskerSignalTrials\whiskerSignalTrial_JF4004_091707_S0_trial90.mat','ws')
load('C:\dan\analysis\behav\analysis\behaviorPerSe\whiskerSignalTrials\whiskerSignalTrial_JF4004_091707_S0_trial90.mat','ws')



%% WhiskerTrialArray and WhiskerSignalTrialArray examples:

%---- DO96 ----------------------------------------------------------------
fnbase = 'K:\LCA\JF8632_070708_DO96\JF8632_070708_DO96_'; mouse_name = 'JF8632'; session_name = 'JF8632_070708_DO96';
trajectory_nums = {[0 1 2],{'D4','D3','D2'}}; % D4=0,D3=1,D2=2
load('C:\dan\analysis\behav\analysis\cellAttached\tifFileNumToBehavTrialNum\DO96.mat','tifFileNums','behavTrialNums','trialTypes') %
w = Whisker.WhiskerTrialArray(fnbase, tifFileNums, behavTrialNums, trajectory_nums, mouse_name, session_name, trialTypes);
% save('C:\dan\analysis\behav\analysis\cellAttached\whiskerTrialArrays\whiskerTrialArray_DO96.mat','w','-v7') 
% load('C:\dan\analysis\behav\analysis\cellAttached\whiskerTrialArrays\whiskerTrialArray_DO96.mat','w') 

xroiPos0=[120 130]; xroiCurv0=[85 130]; 
xroiPos1=[120 130]; xroiCurv1=[85 130]; 
xroiPos2=[120 130]; xroiCurv2=[85 130]; 
xroi = {trajectory_nums{1},[xroiPos0 xroiCurv0], [xroiPos1 xroiCurv1], [xroiPos2 xroiCurv2]};
ws = Whisker.WhiskerSignalTrialArray(w,xroi);

% xroi=[75 131];
% ws = Whisker.WhiskerSignalTrialArray(w,xroi);
% save('C:\dan\analysis\behav\analysis\cellAttached\whiskerSignalTrialArrays\whiskerSignalTrialArray_DO96.mat','ws') 
load('C:\dan\analysis\behav\analysis\cellAttached\whiskerSignalTrialArrays\whiskerSignalTrialArray_DO96.mat','ws') 

%% Plot time projection of whisker for choosing X-pixel ROI:
x = w.trials{28};
figure('Color','white');
x.plot_whisker_time_projection(2, 'b-')
x.plot_whisker_time_projection(1, 'r-')
x.plot_whisker_time_projection(0, 'g-')

%% Plot whisker position:

% Load ws, a WhiskerSignalTrial, then:
[y,t] = ws.get_position(0); 
plot(t,y, 'b-','LineWidth',lw)





