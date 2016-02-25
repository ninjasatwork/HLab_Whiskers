function [trialNums, spikeTimes] = plot_traces_raw_spikes_lfp(cellnum, code, sweepnums, thresh, ...
                FFTWindowSize, FFTStepSize, FFTFreqMin, FFTFreqMax, ylm)

sweepTimes = [];
spikeRate = [];
LFPSpectrogram = cell(1,length(sweepnums));
spikeTimes =  cell(1,length(sweepnums));
trialNums =  zeros(1,length(sweepnums));

sampleRate=10000;
HighPassCutOffInHz = 500;

bitTime = 2; gapTime = 5; preTime = 10; % all in ms, for reading bit code

Wn = HighPassCutOffInHz / (sampleRate/2);
[b,a]=butter(2,Wn,'high');
% ymax = 10; ymin=-10;

if length(thresh)==1
    thresh = repmat(thresh, length(sweepnums), 1);
end

figure('Color','white')
subplot(311)

fs=14;

n=1;
for k=sweepnums
    
    clf   
    [timeInMin, x, bitCode] = CA_load_trace(cellnum,code,k,1:2);

    %-- Plot raw signal:
    subplot(411)
    plot((1:length(x))/sampleRate, x,'k-'); 
    set(gca,'TickDir','out','Box','off','FontSize',fs)
    ylabel('mV','FontSize',fs)
    
    %-- Plot high-pass filtered signal and mark detected spikes:
    subplot(412)
    xx = filtfilt(b,a,x);

    % Use threshold from Quian Quiroga et al. (2004):
    % NEED TO MAKE THIS MORE LOCAL TO AVOID MOVEMENT ARTIFACTS!
%     thresh = 5*median(abs(xx)./0.6745); st = cell_attached_find_spikes(-xx(1:(3*sampleRate)), thresh); % Take just the first 3 s.
%     thresh = 5*median(abs(xx)./0.6745); st = cell_attached_find_spikes(abs(xx), thresh);
%         st = cell_attached_find_spikes(abs(xx), thresh(n));
    st = cell_attached_find_spikes(xx, thresh(n));
%     st = cell_attached_find_spikes(xx(1:(3*sampleRate)), thresh(n)); % Take just the first 3 s for now.
%     st = cell_attached_find_spikes(-xx(1:(3*sampleRate)), thresh(n)); % Take just the first 3 s for now.
        
    plot((1:length(xx))/sampleRate, xx,'k-'); hold on; 
    ylabel('mV','FontSize',fs)
    
    if isempty(ylm)
        ylm = get(gca,'YLim'); ymax = ylm(2);
    else
        ymin = ylm(1); ymax = ylm(2); ylim([ymin ymax]); 
    end
    if st>0
        plot(st/sampleRate,repmat(ymax-.1*ymax, size(st)),'r*')
    end
    set(gca,'TickDir','out','Box','off','FontSize',fs)
    
    %-- Plot sliding-window FFT:
    subplot(413)
    [f,t,YMat] = LCA.FFT_sliding_window(x, sampleRate, FFTWindowSize, FFTStepSize);
    [f,t,YMat] = CA_plot_FFT_sliding_window(f,t,YMat,FFTFreqMin,FFTFreqMax,FFTWindowSize);
    set(gca,'TickDir','out','Box','off','FontSize',fs)
    ylabel('Hz','FontSize',fs)
    xlabel('Sec','FontSize',fs)
    %----------------------------------
    
    subplot(414)
    plot((1:length(bitCode))/sampleRate, bitCode, 'k-')
    
    
    
    if st==0
        FR = 0;
    else
        FR = length(st)/(length(xx)/sampleRate);
    end
%     spikeRate = [spikeRate; FR];
%     sweepTimes = [sweepTimes; timeInMin];
%     LFPSpectrogram{n} = YMat; 
%     
    subplot(411)
%     title([cellnum code ', Sweep=' int2str(k) ', ' num2str(FR) 'Hz, thresh=' num2str(thresh(n))]);
    
    trialNum = read_bit_code(bitCode, sampleRate, bitTime, gapTime, preTime);
    title([cellnum code ', Sweep=' int2str(k) ', Solo trial num=' num2str(trialNum) ',' num2str(FR) 'Hz']);
    
    if st>0
        spikeTimes{n} = st;
    end
    trialNums(n) = trialNum;
    
    n=n+1;
    
    pause
end
% figure; plot(sweepTimes-sweepTimes(1),spikeRate, 'k*')
% xlabel('Min');
% % xlabel('Trial');
% ylabel('Hz')
