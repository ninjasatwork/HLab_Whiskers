function spikeWaveforms = CA_thresh_detect_all_multiunit_spikes(cellnum, code, sweepnums,sampleRate,waveFormTimeInMs, thresh)
%
%
%
%
%
%
% DHO, 4/08.
%



sweepTimes = [];
spikeRate = [];
spikeWaveforms = cell(length(sweepnums),3); % Each entry: {sweepBeginingTime (1 x 1), spikeTimes (1 x n), spikeWaveforms (n x waveFormSamples)}


HighPassCutOffInHz = 500;
Wn = HighPassCutOffInHz / (sampleRate/2);
% [b,a]=butter(2,Wn,'high');
[b,a]=butter(3,Wn,'high');

waveFormSamples = waveFormTimeInMs * sampleRate/1000;


if length(thresh)==1
    thresh = repmat(thresh, length(sweepnums), 1);
end


figure('Color','white')
subplot(311)

fs=14;

n=1;
for k=sweepnums

    clf
    [timeInMin, x, bitCode] = CA_load_trace(cellnum,code,k,1:2); % timeInMin is *NOT* updated for external triggering

    %-- Plot raw signal:
    subplot(311)
    plot((1:length(x))/sampleRate, x,'k-');
    set(gca,'TickDir','out','Box','off','FontSize',fs)
    ylabel('mV','FontSize',fs)

    %-- Plot high-pass filtered signal and mark detected spikes:
    subplot(312)
    xx = filtfilt(b,a,x);

    % Use threshold from Quian Quiroga et al. (2004):
    % NEED TO MAKE THIS MORE LOCAL TO AVOID MOVEMENT ARTIFACTS!
%     thresh = 5*median(abs(xx)./0.6745); st = cell_attached_find_spikes(-xx(1:(3*sampleRate)), thresh); % Take just the first 3 s for now.

%     st = cell_attached_find_spikes(xx, thresh(n));
    st = cell_attached_find_spikes(abs(xx), thresh(n));
%     st = cell_attached_find_spikes(xx(1:(3*sampleRate)), thresh(n)); % Take just the first 3 s for now.
%     st = cell_attached_find_spikes(-xx(1:(3*sampleRate)), thresh(n)); % Take just the first 3 s for now.
        


    if st==0
        wv = [];
    else
        nSpikes = length(st);
        wv = zeros(nSpikes,waveFormSamples);
        for j=1:nSpikes
            if mod(waveFormSamples,2)==0 %even number
                beginSamp = st(j) - (waveFormSamples/2 - 1);
                endSamp = st(j) + waveFormSamples/2;
            else % odd number
                beginSamp = st(j) - (waveFormSamples-1)/2;
                endSamp = st(j) + (waveFormSamples-1)/2;
            end
            if beginSamp > 0 && endSamp <= length(xx)
                wv(j,:) = xx(beginSamp:endSamp);
            else
                disp('Skipping spike at edge of trace bc not enough waveform samples')
            end
            %             plot(wv(j,:))
        end
    end

    spikeWaveforms{k,1} = (k-1)*length(x)/sampleRate; % for now, ignore gaps between trials in setting time %timeInMin;
    spikeWaveforms{k,2} = st./sampleRate; % time in sec
    spikeWaveforms{k,3} = wv;



    plot((1:length(xx))/sampleRate, xx,'k-'); hold on;
    ylm = get(gca,'YLim'); ymax = ylm(2);
    ylabel('mV','FontSize',fs)
    if st>0
        plot(st/sampleRate,repmat(ymax-.1*ymax, size(st)),'r*')
    end
    set(gca,'TickDir','out','Box','off','FontSize',fs)

    
    subplot(313)
    plot((1:length(bitCode))/sampleRate, bitCode, 'k-')
    
    
    %----------------------------------

    %     if st==0
    %         FR = 0;
    %     else
    %         FR = length(st)/(length(xx)/sampleRate);
    %     end
    %     spikeRate = [spikeRate; FR];
    %     sweepTimes = [sweepTimes; timeInMin];
    %     LFPSpectrogram{n} = YMat;
    %
    subplot(311)
    title([cellnum code ', Sweep=' int2str(k)]);
    pause

    n=n+1;
end

