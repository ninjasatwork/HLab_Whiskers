d = 'C:\Documents and Settings\oconnord\My Documents\MATLAB\testOneByOneConstruction';
trajectory_nums = [0 1];

% Example use to process all .whiskers files in directory d:
Whisker.makeAllDirectory_WhiskerTrial(d,trajectory_nums,'faceSideInImage','top',...
    'protractionDirection','rightward','barRadius',3,'pxPerMm',13.74,'imagePixelDimsXY',[180 250],'mouseName','JF25399',...
    'sessionName','091101-JF25399b');

Whisker.makeAllDirectory_WhiskerSignalTrial(d,'polyRoiInPix',[25 75],'follicleExtrapDistInPix',13); % If want to calculate forces, use ''follicleExtrapDistInPix' argument.
% Whisker.makeAllDirectory_WhiskerSignalTrial(d,'polyRoiInPix',[0 200]); % Otherwise don't; it will be faster.

Whisker.makeAllDirectory_WhiskerTrialLite(d,'r_in_mm',3.6);
Whisker.makeAllDirectory_WhiskerTrialLite(d,'r_in_mm',3.6,'calc_forces',true);

wl = Whisker.WhiskerTrialLiteArray(d);
% save('test_WTLA.mat','wl');

% View the resulting timeseries:
tid=0; contact_tid = 1;
Whisker.view_WhiskerTrialLiteArray(wl,tid,contact_tid)
% Or, without contact annotations:
Whisker.view_WhiskerTrialLiteArray(wl,tid)

%%

tidPairs = [10 0; 11 1; 12 2]; %An Nx2 matrix where each row is of form: [tid_to_copy tid_target]
Whisker.allDirectory_copyTraj(d,'include_files',{'LTPJF41501_081909_C1_B_0001'},'trajectoryIDPairs',tidPairs,'follicleExtrapDistInPix',13);

% Whisker.makeAllDirectory_WhiskerTrialLite(d,'include_files',{'LTPJF41501_081909_C1_B_0001'},'r_in_mm',3.6);
Whisker.makeAllDirectory_WhiskerTrialLite(d,'include_files',{'LTPJF41501_081909_C1_B_0001'},'r_in_mm',3.6,'calc_forces',true);

%%
tid=0;
figure
subplot(211)
wl.plot_whisker_angle(tid,'k-')
subplot(212)
wl.plot_whisker_curvature(tid,'k-')

%%
% Load a WhiskerSignalTrial object named ws, then:
tid = 0; r_in_mm = 1; near_base_in_mm = 0.2;
[thetap,kappap,y,x,t] = ws.get_theta_kappa_at_roi_point(tid,r_in_mm);
[thetaNearBase,kappaNearBase,y,x,t] = ws.get_theta_kappa_at_roi_point(tid,near_base_in_mm);

figure 
subplot(211)
plot(t,thetaNearBase,'k-'); xlabel('Sec'); ylabel('Theta at base(deg)')
subplot(212); plot(t,kappap,'k-'); xlabel('Sec'); ylabel('Kappa (1/pixels)')


%%
wl = Whisker.WhiskerTrialLite(ws,'r_in_mm',2);


%% Alternative uses:
% Whisker.makeAllDirectory_WhiskerSignalTrial(d,'polyRoiInPix',[0 200],'follicleExtrapDistInPix',13);

% Example use to process only those .whiskers files in directory d specified in argument
% 'include_files':
Whisker.makeAllDirectory_WhiskerTrial(d,trajectory_nums,'include_files',{'091101-JF25399b_005',...
    '091101-JF25399b_006','091101-JF25399b_008'},...
    'faceSideInImage','top',...
    'barRadius',3,'pxPerMm',13.74,'imagePixelDimsXY',[180 250],'mouseName','JF25399',...
    'sessionName','091101-JF25399b');

Whisker.makeAllDirectory_WhiskerTrial(d,trajectory_nums,'include_files',{'091101-JF25399b_005',...
    '091101-JF25399b_008','091101-JF25399b_007','091101-JF25399b_006'},...
    'faceSideInImage','top',...
    'barRadius',3,'pxPerMm',13.74,'imagePixelDimsXY',[180 250],'mouseName','JF25399',...
    'sessionName','091101-JF25399b');

Whisker.makeAllDirectory_WhiskerSignalTrial(d,'include_files',{'091101-JF25399b_005',...
    '091101-JF25399b_008','091101-JF25399b_007','091101-JF25399b_006'},'polyRoiInPix',[0 200]);


%%  Example use for Leo for getting (x,y) at tip:
tid = 0;
[theta,kappa,y,x,t] = ws.get_theta_kappa_at_point(tid,'max');


%% Plot whiskers from all frames, radial limits colored:

tid = 0;
figure
ws.plot_fitted_whisker_time_projection(tid, 'k-'); hold on
ws.plot_fitted_whisker_ROI_time_projection(tid, 'r-')
% set(gca,'XLim',[0 360],'YLim',[0 320])
set(gca,'XLim',[0 180],'YLim',[0 250])













