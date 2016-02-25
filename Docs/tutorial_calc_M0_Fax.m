function g = tutorial_calc_M0_Fax(ws,tid,contact_tid,r_point,varargin)
%
%
% INPUTS:
%
%   ws: A WhiskerSignalTrialArray.
%   tid: The trajectoryID of the whisker to make measurements on.
%   contact_tid: The trajectory ID of the whisker used to score contacts.
%   r_point: Radial distance along whisker at which to measure kappa. In mm.
%
% OUTPUTS:
%
%   g: A structure that can be passed to Whisker.view_M0_Fax for viewing.
%      Has the following fields:
%           g.frames
%           g.M0
%           g.Faxial
%           g.dkappa
%           g.distanceToPoleCenter
%           g.thetaAtContact
%           g.thetaAtBase
%           g.contactFrames
%           g.fileNames
%           g.follX
%           g.follY
%
%
%
whisker_radius_at_base = 33.5; % In microns. 
whisker_length = 16; % In mm.
youngs_modulus = 5e9; % In Pa
baseline_time_end = 0.1; % In sec.


[M0,Faxial,t,dkappa,Fnorm,thetaAtBase,thetaAtContact,distanceToPole] = ...
    ws.calc_M0_Faxial(tid,r_point,whisker_radius_at_base, whisker_length,youngs_modulus,baseline_time_end);

[follY,follX,jnk] = ws.get_cached_follicle_coords(tid);

contactFrames = cellfun(@(x) x.time{x.trajectoryIDs==contact_tid}/x.framePeriodInSec, ws.trials,'UniformOutput',false);
fileNames = cellfun(@(x) x.trackerFileName, ws.trials,'UniformOutput',false);

framePeriodInSec = ws.trials{1}.framePeriodInSec; % Assume frame rate is same for all trials.
frames = cellfun(@(x) x/framePeriodInSec, t,'UniformOutput',false); % Convert time to frames by dividing by frame period of 0.002 sec.

g.frames = frames;
g.M0 = M0;
g.Faxial = Faxial;
g.dkappa = dkappa;
g.distanceToPole = distanceToPole;
g.thetaAtContact = thetaAtContact;
g.thetaAtBase = thetaAtBase;
g.contactFrames = contactFrames;
g.fileNames = fileNames;
g.follX = follX;
g.follY = follY;

















