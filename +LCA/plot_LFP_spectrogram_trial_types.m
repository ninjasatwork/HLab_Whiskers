t = s.trialNums;

h = intersect(t, b.hitTrialNums);
m = intersect(t, b.missTrialNums);
fa = intersect(t, b.falseAlarmTrialNums);
cr = intersect(t, b.correctRejectionTrialNums);

[S_h, count_h] = s.get_LFP_spectrogram_mean(h);
[S_cr, count_cr] = s.get_LFP_spectrogram_mean(cr);
[S_fa, count_fa] = s.get_LFP_spectrogram_mean(fa);
[S_m, count_m] = s.get_LFP_spectrogram_mean(m);

freqMin=0; freqMax=30; FFTWindowSize=1; xmin=0; xmax=5;
figure('Color','white')
spct = S_h; subplot(221)
if ~isempty(S_h)
    LCA.plot_FFT_sliding_window(spct.f, spct.t, spct.S, freqMin, freqMax, FFTWindowSize); 
    xlim([xmin xmax]); title(['H, n=' int2str(count_h)])
%     colorbar
    clim_h = get(gca,'CLim');
else
    clim_h = NaN;    
end

spct = S_cr; subplot(222)
if ~isempty(S_cr)
    LCA.plot_FFT_sliding_window(spct.f, spct.t, spct.S, freqMin, freqMax, FFTWindowSize); 
    xlim([xmin xmax]); title(['CR, n=' int2str(count_cr)])
%     colorbar
    clim_cr = get(gca,'CLim');
else
    clim_cr = NaN;
end

if ~isempty(S_fa)
    spct = S_fa; subplot(223) 
    LCA.plot_FFT_sliding_window(spct.f, spct.t, spct.S, freqMin, freqMax, FFTWindowSize); 
    xlim([xmin xmax]); title(['FA, n=' int2str(count_fa)])
%     colorbar
    clim_fa = get(gca,'CLim');
else
    clim_fa = NaN;
end

if ~isempty(S_m)
    spct = S_m; subplot(224) 
    LCA.plot_FFT_sliding_window(spct.f, spct.t, spct.S, freqMin, freqMax, FFTWindowSize); 
    xlim([xmin xmax]); title(['M, n=' int2str(count_m)])
%     colorbar
    clim_m = get(gca,'CLim');
else
    clim_m = NaN;
end

% Set all to have same color scale:
% subplot(221); clim_h = get(gca,'CLim');
% subplot(222); clim_cr = get(gca,'CLim');
% subplot(223); clim_fa = get(gca,'CLim');
% subplot(224); clim_m = get(gca,'CLim');

% clim_max = nanmax(nanmax([clim_h; clim_cr; clim_fa; clim_m]));
% clim_min = nanmin(nanmin([clim_h; clim_cr; clim_fa; clim_m]));

clim_max = nanmax(nanmax([clim_h; clim_cr;])); % too few FAs and Misses sometimes.
clim_min = nanmin(nanmin([clim_h; clim_cr;])); % too few FAs and Misses sometimes.

% clim_max = nanmax(nanmax([clim_cr;])); 
% clim_min = nanmin(nanmin([clim_cr;])); % 

subplot(221); set(gca,'CLim', [clim_min clim_max]);
subplot(222); set(gca,'CLim', [clim_min clim_max]);
subplot(223); set(gca,'CLim', [clim_min clim_max]);
subplot(224); set(gca,'CLim', [clim_min clim_max]);









