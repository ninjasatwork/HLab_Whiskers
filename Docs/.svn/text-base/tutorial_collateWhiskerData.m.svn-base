
%% Build the WhiskerSignalTrialArray from raw .whiskers files:
%

trajectoryIDs = [0 1];
fn = {'091101-JF25399b_005', '091101-JF25399b_006', '091101-JF25399b_007', '091101-JF25399b_008',trajectoryIDs};
w = Whisker.WhiskerTrialArray(fn); w.set_imagePixelDimsXY([180 250]); w.set_faceSideInImage('top'); 
w.set_pxPerMm(17.16); w.fit_polys;  w.set_faceSideInImage('top');
% save('whiskerTrialArray_091101-JF25399b_entireTrial_small','w')
% load('whiskerTrialArray_091101-JF25399b_entireTrial_small','w')

roiPos0=[1 18]; roiCurv0=[1 40]; extrap_distance_in_pix = 13;
ws = Whisker.WhiskerSignalTrialArray(w,[roiPos0 roiCurv0]); % WhiskerSignalTrialArray comprises WhiskerSignalTrial
ws.recompute_cached_follicle_coords(extrap_distance_in_pix,0);

% save('whiskerSignalTrialArray_091101-JF25399b_entireTrial_small','ws')
% load('whiskerSignalTrialArray_091101-JF25399b_entireTrial_small','ws')



%% Measurements inside a separate polynomial fit over constant region of arc-length:

trajectoryIDs = [0 1];
fn = {'091101-JF25399b_005', '091101-JF25399b_006', '091101-JF25399b_007', '091101-JF25399b_008',trajectoryIDs};
w = Whisker.WhiskerTrialArray(fn); w.set_imagePixelDimsXY([180 250]); w.set_faceSideInImage('top'); w.set_pxPerMm(17.16); 

% [xx,yy] = ginput; plot(xx,yy,'k-'); disp(xx'); disp(yy')
xx = [43.4879   60.1815   79.7782   93.0847  105.1815  116.6734];
yy = [34.7222   45.0708   45.3431   44.5261   43.1645   34.9946];
w = w.set_mask_from_points([0 1],xx,yy);

w.fit_polys_roi([0 100]); % In pixels 
% save('whiskerTrialArray_091101-JF25399b_entireTrial_small','w')
% load('whiskerTrialArray_091101-JF25399b_entireTrial_small','w')

roiPos0=[1 5]; roiCurv0=[1 40]; extrap_distance_in_pix = 13;
ws = Whisker.WhiskerSignalTrialArray(w,[roiPos0 roiCurv0]); 
ws.recompute_cached_follicle_coords(extrap_distance_in_pix,0);

figure; k=1;
ws.trials{k}.plot_fitted_whisker_time_projection(0, 'k-'); hold on
ws.trials{k}.plot_fitted_whisker_ROI_time_projection(0, 'r-')

% save('whiskerSignalTrialArray_091101-JF25399b_entireTrial_small','ws') 
load('whiskerSignalTrialArray_091101-JF25399b_entireTrial_small','ws')

%---
tid=0; contact_tid = 1;
r_point = 2.33;  % In mm.

g = tutorial_calc_M0_Fax(ws,tid,contact_tid,r_point);
% save('M0_Faxial_091101-JF25399b-2.mat','g')
% load('M0_Faxial_091101-JF25399b-2')

Whisker.view_M0_Fax(g)





%% Example plotting M0, Faxial, etc:
%
load('whiskerSignalTrialArray_091101-JF25399b_entireTrial_small','ws')
tid=0; contact_tid = 1;
r_point = 2.33;  % In mm.

g = tutorial_calc_M0_Fax(ws,tid,contact_tid,r_point);
% save('M0_Faxial_091101-JF25399b.mat','g')
load('M0_Faxial_091101-JF25399b')

Whisker.view_M0_Fax(g)


%% - Test calc_M0_Faxial():
load('whiskerSignalTrialArray_091101-JF25399b_entireTrial_small','ws')
tid=0; contactTid = 1;
x = ws.trials{1}; 
x = ws.trials{2}; 

contactFrames = x.time{x.trajectoryIDs==contactTid}/x.framePeriodInSec;


r_point = 2.33; %3.6; % In mm.
whisker_radius_at_base = 33.5; % In microns. 
whisker_length = 16; % In mm.
youngs_modulus = 5e9; % In Pa
baseline_time_end = 0.1; % In sec.

[M0,Faxial,t,dkappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPoleCenter] = ...
    x.calc_M0_Faxial(tid,r_point,whisker_radius_at_base, whisker_length,youngs_modulus,baseline_time_end);

f = t/x.framePeriodInSec;

figure('Position',[-1286 22 1098 1081]);
subplot(611); plot(f,distanceToPoleCenter,'r.-'); title('Distance to pole center'); ylabel('mm')
subplot(612); plot(f,thetaAtBase,'k.-'); hold on; plot(f,thetaAtContact,'r.-'); title('Theta at base (red) and contact (red)'); ylabel('Deg')
subplot(613); plot(f,dkappa,'k.-'); title(['Kappa at ' num2str(r_point) ' mm']); ylabel('1/mm')
subplot(614); plot(f,Fnorm,'k.-'); title('|F|'); ylabel('N')
subplot(615); plot(f,Faxial,'k.-'); title('F_{axial}'); ylabel('N')
subplot(616); plot(f,M0,'k.-'); title('M_0'); ylabel('Newton-meters')
for k=1:6
    subplot(6,1,k)
    xlim([0 749])
    ylm = get(gca,'YLim');
    for q=1:length(contactFrames)
        line([contactFrames(q) contactFrames(q)],ylm,'Color','b')
    end
end
          


%% Example of building Whisker(Signal)TrialArray from raw data with 3 tracked whiskers:
%
% fnbase = 'K:\LCA\JF8632_062308_DO79\JF8632_062308_DO79_'; mouse_name = 'JF8632'; session_name = 'JF8632_062308_DO79';
trajectory_nums = {[0 1 2],{'D4','D3','D2'}};  % D4=0,D3=1,D2=2
% load('C:\dan\analysis\behav\analysis\cellAttached\tifFileNumToBehavTrialNum\DO79.mat','tifFileNums','behavTrialNums','trialTypes') % 51:157, 58:164
w = Whisker.WhiskerTrialArray(fnbase, tifFileNums, behavTrialNums, trajectory_nums, mouse_name, session_name, trialTypes);
% save('C:\dan\analysis\behav\analysis\cellAttached\whiskerTrialArrays\whiskerTrialArray_DO79.mat','w','-v7') 
% load('C:\dan\analysis\behav\analysis\cellAttached\whiskerTrialArrays\whiskerTrialArray_DO79.mat','w') 

roiPos0=[1 20]; roiCurv0=[1 40]; 
roiPos1=[1 12]; roiCurv1=[1 40]; 
roiPos2=[1 12]; roiCurv2=[1 40]; 
roi = {trajectory_nums{1},[roiPos0 roiCurv0], [roiPos1 roiCurv1], [roiPos2 roiCurv2]};

% ws = Whisker.WhiskerSignalTrialArray(w,roi);

% [xx,yy] = ginput; plot(xx,yy,'k-'); disp(xx'); disp(yy')
xx = [48.3180   58.2719   71.1290   84.4009  101.4055  125.0461];
yy = [23.0263   28.1433   28.1433   27.4123   27.4123    9.1374];

ws = ws.set_mask_from_points([0 1],xx,yy);
ws = ws.recompute_cached_mean_theta_kappa(roi);

ws = ws.set_mask_from_points(trajectory_nums{1},xx,yy);
ws = ws.recompute_cached_mean_theta_kappa(roi);


%% Plot whiskers from all frames, radial limits colored:

x = ws.trials{1}; 
x = ws.trials{25}; 
x = ws.trials{30}; 
x = ws.trials{45}; 

figure('Color','white'); hold on
% x.plot_fitted_whisker_time_projection(2, 'b-',[],{roiPos2,'k-'})
% x.plot_fitted_whisker_time_projection(1, 'r-',[],{roiPos1,'k-'})
x.plot_fitted_whisker_time_projection(0, 'g-',[],{roiPos0,'k-'})
% x.plot_fitted_whisker_time_projection(3, 'm-',[],{roiPos3,'k-'})
set(gca,'XLim',[0 180],'YLim',[0 250])
x.plot_mask(0,'k-');

figure('Color','white'); hold on
x.plot_fitted_whisker_time_projection(2, 'b-',[],{roiCurv2,'k-'})
x.plot_fitted_whisker_time_projection(1, 'r-',[],{roiCurv1,'k-'})
x.plot_fitted_whisker_time_projection(0, 'g-',[],{roiCurv0,'k-'})
% x.plot_fitted_whisker_time_projection(3, 'm-',[],{roiCurv3,'k-'})
set(gca,'XLim',[0 180],'YLim',[0 250])
x.plot_mask(0,'k-');

%% Plot theta and kappa:
x = ws.trials{1}; 

figure
subplot(211); hold on; 
x.plot_whisker_angle(0,'b-'); 
% x.plot_whisker_angle(1,'g-');
% x.plot_whisker_angle(2,'m-');
% x.plot_whisker_angle(3,'k-');
ylabel('\theta (deg)')
subplot(212); hold on
x.plot_whisker_curvature(0,'b-'); 
% x.plot_whisker_curvature(1,'g-'); 
% x.plot_whisker_curvature(2,'m-'); 
% x.plot_whisker_curvature(3,'k-'); 
ylabel('\kappa (1/pixels)')








