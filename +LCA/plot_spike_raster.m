function plot_spike_raster(spike_times, option_string)
%
% plot_spike_raster(spike_times, 'option_string')
%
%   spike_times:  Object of type LCA.SpikesTrialArray.
%   option_string: Either 'BehavTrialNum' or 'Sequential'.
%                  'BehavTrialNum' will plot behavioral trial number on
%                   y-axis (including any gaps from missing trials).
%                   'Sequential' will just plot 1:length(number of trials)
%                   on the y-axis.
%
% DHO, 5/08.
%

if ~isa(spike_times, 'LCA.SpikesTrialArray')
    error('First argument must be of type LCA.SpikesTrialArray.')
end

hold on

allSpikeTimes = [];
switch option_string

    case 'BehavTrialNum'

        for k=1:length(spike_times)
            t = spike_times.spikesTrials{k}.spikeTimes;
            if t > 0
%                 plot(t./spike_times.spikesTrials{k}.sampleRate, spike_times.spikesTrials{k}.trialNum, 'k.')
                trialSpikeTimes = t./spike_times.spikesTrials{k}.sampleRate;
                allSpikeTimes = [allSpikeTimes; repmat(spike_times.spikesTrials{k}.trialNum, ...
                                    size(trialSpikeTimes)), trialSpikeTimes];
            end
        end
        plot(allSpikeTimes(:,2), allSpikeTimes(:,1), 'k.')
        ylabel('Behavior trial number','FontSize',15)
        xlabel('Sec','FontSize',15)
        
    case 'Sequential'

        for k=1:length(spike_times)
            t = spike_times.spikesTrials{k}.spikeTimes;
            if t > 0
%                 plot(t./spike_times.spikesTrials{k}.sampleRate, k, 'k.')
                trialSpikeTimes = t./spike_times.spikesTrials{k}.sampleRate;
                allSpikeTimes = [allSpikeTimes; repmat(k, size(trialSpikeTimes)), trialSpikeTimes];
            end
        end
        plot(allSpikeTimes(:,2), allSpikeTimes(:,1), 'k.')
        ylabel('Sequential trial number','FontSize',15)
        xlabel('Sec','FontSize',15)

    otherwise
        error('Invalid string argument')

end