t = T.trialNums; 

h = T.hitTrialNums;
m = T.missTrialNums;
fa = T.falseAlarmTrialNums;
cr = T.correctRejectionTrialNums;


figure('Color','white'); fs = 10;

subplot(221)
T.plot_spike_raster(h,'Sequential')
title(['H, n=' int2str(length(h))])
% set(gca,'FontSize',fs)
ylabel('Sequential trial number','FontSize',10)
xlabel('Sec','FontSize',10)
xlim([0 5])

subplot(222)
T.plot_spike_raster(cr,'Sequential')
title(['CR, n=' int2str(length(cr))])
% set(gca,'FontSize',fs)
ylabel('Sequential trial number','FontSize',10)
xlabel('Sec','FontSize',10)
xlim([0 5])

subplot(223)
T.plot_spike_raster(fa,'Sequential')
title(['FA, n=' int2str(length(fa))])
% set(gca,'FontSize',fs)
ylabel('Sequential trial number','FontSize',10)
xlabel('Sec','FontSize',10)
xlim([0 5])

subplot(224)
if ~isempty(m)
    T.plot_spike_raster(m,'Sequential')
end
title(['M, n=' int2str(length(m))])
% set(gca,'FontSize',fs)
ylabel('Sequential trial number','FontSize',10)
xlabel('Sec','FontSize',10)
xlim([0 5])





