%
% Transform an XSG file into a more useful data structure.
%
% Includes dependent properties for filtered versions of signal to save
% disk space.
%
% Have "primary" cell that's patched.  Have "secondary" voltage trace after
% subtracting off spike waveforms from primary cell.  Can't use single
% (low) threshold for both because each primary spike has multiple peaks
% that will cross threshold. 
%
%
%
% DHO, 4/08.
%
classdef Sweep < handle

    properties

        cellNum = '';
        cellCode = '';
        trialNum  = []; % as parsed from the bitCode.
        xsgFileNum = [];
        
        bitCode   = []; % truncated from bitCode XSG channel

        useFlag   = 1; % set to 0 to mark bad sweeps that should be trimmed.
        timeStamp = [];
        sampleRate = []; %10000;
        rawSignal = [];

        bitCodeBitTime = 2; % ms
        bitCodeGapTime = 5; % ms
        bitCodePreTime = 10; % ms, Should make these as optional arguments to constructor.

        highPassCutOffInHz = 600;
        LFPBandPassCutOffsInHz = [2 100];

        FFTWindowSize = 1; % sec, for LFP spectrogram.
        FFTStepSize = 0.1; % sec, for LFP spectrogram.
        LFPSpectrogramDisplayFreqMin = 2; % Hz, for display only.
        LFPSpectrogramDisplayFreqMax = 30; % Hz, for display only.

        spikeThreshold = 1.5; % mV % Will need to sometimes use abs() sometimes not.
        spikeThresholdSecondary = 0.2; % mV 
        artifactThreshold = -5 % mV 
        
        waveformTimeInMs = 3;%6;
        
    end

    properties (Dependent = true)
        highPassFilteredSignal
        spikeTimes
        spikeWaveforms

        sweepLengthInSamples
        
        LFP
        LFPSpectrogram

        highPassFilteredSignalSecondary % after subtraction of primary spikes.
        spikeTimesSecondary
        spikeWaveformsSecondary
        
    end


    % Needed methods: loading/saving, exporting concatentated (primary and secondary) waveforms for input
    % to MClust.
    methods (Access = public)
        function obj = Sweep(cellnum,cellcode,sweepnum)
            if nargin > 0
                [timeInMin, obj.sampleRate, obj.rawSignal, obj.bitCode] = LCA.load_xsg_trace(cellnum, cellcode, sweepnum, 1:2);
%                 if length(obj.bitCode) > (0.08*obj.sampleRate) %800
%                     obj.bitCode = obj.bitCode(1:(0.08*obj.sampleRate)); %code(1:800);
%                 end
                obj.trialNum = LCA.read_bit_code(obj.bitCode, obj.sampleRate, ...
                    obj.bitCodeBitTime, obj.bitCodeGapTime, obj.bitCodePreTime);
                obj.cellNum = cellnum;
                obj.cellCode = cellcode;
                obj.xsgFileNum = sweepnum;
            end
        end

%         function delete(obj)
%             % Don't need any special cleaning up.
%         end
        
%         function disp(obj)
% %             sprintf('bitCode: %g\nSample Number: %g\nModulus: %1.5g\n',...
% %                 obj.bitCode, obj.trialNum, obj.useFlag, obj.timeStamp,...
% %                 obj.sampleRate, size(obj.rawSignal));
%                 sprintf(['Object LCA.Sweep:\n\nbitCode\ntrialNum\nuseFlag\n' ...
%                     'timeStamp\nsampleRate\nrawSignal\nbitCodeBitTime\n'])
%         end
%         
        function plot_LFP_spectrogram(obj)
            LCA.plot_FFT_sliding_window(obj.LFPSpectrogram.f, ...
                obj.LFPSpectrogram.t,obj.LFPSpectrogram.S, ...
                obj.LFPSpectrogramDisplayFreqMin, ...
                obj.LFPSpectrogramDisplayFreqMax, ...
                obj.FFTWindowSize);
        end

        function plot_spike_times(obj, yval)
            st = obj.spikeTimes;
            if nargin < 2
                yval = 3;
            end
            if st > 0
                plot(st,yval*ones(size(st)),'r*')
            end
        end
        
        function plot_spike_times_secondary(obj, yval)
            st = obj.spikeTimesSecondary;
            if nargin < 2
                yval = 3;
            end
            if st > 0
                plot(st,yval*ones(size(st)),'g*')
            end
        end
        
        function set_spike_threshold(obj, thresh)
            obj.spikeThreshold = thresh; 
        end
        
        function set_spike_threshold_secondary(obj, thresh)
            obj.spikeThresholdSecondary = thresh;
        end
        
        function set_waveform_time_in_ms(obj, ms)
            obj.waveformTimeInMs = ms;
        end
        
        
    end
    

    methods (Access = private, Static = true)
        function [beginSamp, endSamp] = get_waveform_sample_indices(spike_time, number_waveform_samples)
            if mod(number_waveform_samples,2)==0 %even number
                beginSamp = spike_time - (number_waveform_samples/2 - 1);
                endSamp = spike_time + number_waveform_samples/2;
            else % odd number
                beginSamp = spike_time - (number_waveform_samples-1)/2;
                endSamp = spike_time + (number_waveform_samples-1)/2;
            end       
            beginSamp = int32(squeeze(beginSamp));
            endSamp = int32(squeeze(endSamp));
        end
    end

    methods % Dependent property methods; cannot have attributes.
        function value = get.highPassFilteredSignal(obj)
            Wn = obj.highPassCutOffInHz / (obj.sampleRate/2);
            [b,a]=butter(2,Wn,'high');
            value = filtfilt(b, a, obj.rawSignal);
        end

        function value = get.LFP(obj)
            W1 = obj.LFPBandPassCutOffsInHz(1) / (obj.sampleRate/2);
            W2 = obj.LFPBandPassCutOffsInHz(2) / (obj.sampleRate/2);
            [b,a]=butter(2,[W1 W2]);
            value = filtfilt(b, a, obj.rawSignal);
        end

        function value = get.LFPSpectrogram(obj)
            [f,t,S] = LCA.FFT_sliding_window(obj.rawSignal,...
                obj.sampleRate, obj.FFTWindowSize, obj.FFTStepSize);
            value.f = f;
            value.t = t;
            value.S = S;
        end

        function value = get.spikeTimes(obj)
            value = LCA.find_spikes(obj.highPassFilteredSignal,obj.spikeThreshold,obj.artifactThreshold);
        end

        function value = get.spikeWaveforms(obj)
            st = obj.spikeTimes;
            if st==0
                wv = [];
                stWithoutEdgeSpikes = [];
            else
                waveformSamples = obj.waveformTimeInMs * obj.sampleRate/1000;
                nSpikes = length(st);
%                 y = obj.highPassFilteredSignal;
                y = obj.rawSignal;

%                 wv = zeros(nSpikes,waveformSamples);
                wv = []; % Don't yet know how many spikes are at edges of sweep and excludes, so can't initials with nSpikes dimensions. 
                stWithoutEdgeSpikes = []; % Have to trim spike times data structure of spikes that get thrown out due to occuring too close to edge of sweep to get waveform.
                for j=1:nSpikes
                    [beginSamp, endSamp] = LCA.Sweep.get_waveform_sample_indices(st(j), waveformSamples);
                    if beginSamp > 0 && endSamp <= length(y)
                        wv = [wv, y(beginSamp:endSamp)];
                        stWithoutEdgeSpikes = [stWithoutEdgeSpikes; st(j)];
                    else
                        disp('Skipping spike at edge of trace bc not enough waveform samples')
                    end
                end
            end

            value = {};
            value{1} = stWithoutEdgeSpikes./obj.sampleRate; % time in sec
            value{2} = wv;
        end
        
        function value = get.sweepLengthInSamples(obj)
            value = length(obj.rawSignal);
        end      

        function value = get.highPassFilteredSignalSecondary(obj)
            thresholdSpikeTimes = obj.spikeTimes;
            if isempty(thresholdSpikeTimes)
                value = obj.highPassFilteredSignal;
            else
                value = obj.highPassFilteredSignal;
                Y = obj.spikeWaveforms;
%                 w = Y{2};
                waveformSpikeTimes = Y{1}.*obj.sampleRate; % Convert back to indices from ms;
                
                numWaveformSpikes = size(waveformSpikeTimes,1);
                numWaveformSamples = obj.waveformTimeInMs * obj.sampleRate/1000;
                for k=1:numWaveformSpikes
                    [beginSamp, endSamp] = LCA.Sweep.get_waveform_sample_indices(waveformSpikeTimes(k), numWaveformSamples);
%                     value(beginSamp:endSamp) = value(beginSamp:endSamp) - w(k,:)'; % Best to set to 0 to avoid overlaps giving negative values, but doing actual subtraction helps with QA.
                      value(beginSamp:endSamp) = 0;
                end
                numThresholdSpikes = size(thresholdSpikeTimes, 1);
                if numThresholdSpikes > numWaveformSpikes
                    disp('Number of threshold detected spikes is greater than number of spike waveforms')
                    % NEED TO ZERO SAMPLES FOR ANY PRIMARY CELL SPIKES THAT AREN'T
                    % RECORDED IN WAVEFORMS BECAUSE OF BEING AT EDGE.  NOW THEY GET
                    % CAUGHT AS SECONDARY SPIKES.
                    % For now just set to 0 number of samples in waveform
                    % both at the beginning and the end of the sweep. Later
                    % should find edge spikes and set only those to 0.
                    value(1:numWaveformSamples) = 0;
                    value(((end-numWaveformSamples)+1):end) = 0;
                end
            end
        end

        function value = get.spikeTimesSecondary(obj)
            value = LCA.find_spikes(obj.highPassFilteredSignalSecondary,obj.spikeThresholdSecondary,obj.artifactThreshold);
        end
        
        function value = get.spikeWaveformsSecondary(obj)
            st = obj.spikeTimesSecondary;
            if st==0
                wv = [];
                stWithoutEdgeSpikes = [];
            else
                waveformSamples = obj.waveformTimeInMs * obj.sampleRate/1000;
                nSpikes = length(st);
                y = obj.highPassFilteredSignalSecondary;
%                 wv = zeros(nSpikes,waveformSamples);
                wv = []; % Don't yet know how many spikes are at edges of sweep and excludes, so can't initials with nSpikes dimensions. 
                stWithoutEdgeSpikes = []; % Have to trim spike times data structure of spikes that get thrown out due to occuring too close to edge of sweep to get waveform.
                for j=1:nSpikes
                    [beginSamp, endSamp] = LCA.Sweep.get_waveform_sample_indices(st(j), waveformSamples);
                    if beginSamp > 0 && endSamp <= length(y)
                        wv = [wv, y(beginSamp:endSamp)];
                        stWithoutEdgeSpikes = [stWithoutEdgeSpikes; st(j)];
                    else
                        disp('Skipping spike at edge of trace bc not enough waveform samples')
                    end
                end
            end

            value = {};
            value{1} = stWithoutEdgeSpikes./obj.sampleRate; % time in sec
            value{2} = wv;
        end
    end

end




