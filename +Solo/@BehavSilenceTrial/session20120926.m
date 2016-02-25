%  obj = BehavTrialArray_JY(x, session_name)
mouseName   = 'ANM169796';
sessionName = '120926';

fn =['C:\Work\Projects\BehavingVm\Data\Solodata\ANM169796\data_@pole_disc_jyinactivationobj_ANM169796_120926a.mat'];
b = Solo.BehavSilenceTrialArray(fn, sessionName); b.trim = [1 1];  b.performanceRegion = [1 300]; 

%fn =['G:\BehavingVm\Data\Solodata\ANM148967\data_@pole_disc_jyobj_ANM148967_120523a'];
%b = Solo.BehavTrialArray_JY(mouseName, sn); b.trim = [22 35];  b.performanceRegion = [1 125]; % spike rate increasing
figure; b.plot_scored_trials
pc = performance(b)
save (['C:\Work\Projects\BehavingVm\Data\Behaviordata\ANM169796\solo_' mouseName '_' sessionName '.mat'], 'b')
