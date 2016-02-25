function viewer_LFPvsWhisker(varargin)
%
% USAGE:
%  viewer_LFPvsWhisker(s, T)
%
%   s, LCA.SweepArray object.
%   T, LCA.TrialArray object.
%
%   These two inputs must be from the same recording.
%
%   Additional arguments in varargin are reserved for internal, recursive
%       use of this function.
%
%

if nargin==2 && isa(varargin{1},'LCA.SweepArray') && isa(varargin{2},'LCA.TrialArray') % Called for first time, with only SweepArray and TrialArray arguments
    s = varargin{1}; T = varargin{2};
    h=figure('Color','white'); ht = uitoolbar(h);
    a = .20:.05:0.95; b(:,:,1) = repmat(a,16,1)'; b(:,:,2) = repmat(a,16,1); b(:,:,3) = repmat(flipdim(a,2),16,1);
    bbutton = uipushtool(ht,'CData',b,'TooltipString','Back');
    fbutton = uipushtool(ht,'CData',b,'TooltipString','Forward','Separator','on');
    set(fbutton,'ClickedCallback',['LCA.' mfilename '(''next'')'])
    set(bbutton,'ClickedCallback',['LCA.' mfilename '(''last'')'])
    
    m=uimenu(h,'Label','Display Type','Separator','on');
    uimenu(m,'Label','Raw voltage','Callback',['LCA.' mfilename '(''rawSignal'',''none'')'])
    uimenu(m,'Label','High-pass filtered','Callback',['LCA.' mfilename '(''highPassFilteredSignal'',''none'')'])
    uimenu(m,'Label','LFP','Callback',['LCA.' mfilename '(''LFP'',''none'')']);
    uimenu(m,'Label','High-pass filtered + spike times','Callback',['LCA.' mfilename '(''highPassFilteredSignal'',''spikeTimes'')'])
    uimenu(m,'Label','Raw voltage + spike times','Callback',['LCA.' mfilename '(''rawSignal'',''spikeTimes'')'])
    uimenu(m,'Label','LFP + spike times','Callback',['LCA.' mfilename '(''LFP'',''spikeTimes'')'])
    uimenu(m,'Label','High-pass filtered + spike times + threshold','Callback',['LCA.' mfilename '(''highPassFilteredSignal'',''spikeTimesPlusThreshold'')'])
    
    uimenu(m,'Label','Bitcode channel','Callback',['LCA.' mfilename '(''bitCode'',''none'')'])
    
    uimenu(h,'Label','Jump to sweep','Separator','on','Callback',['LCA.' mfilename '(''jumpToSweep'')']);
    
    g = struct('s',s,'T',T,'sweepNum',1,'trialList','','displayType','rawSignal','displayTypeMinor','none');
    
    y=s.sweeps{g.sweepNum}.rawSignal; x=(1:length(y))/s.sweeps{g.sweepNum}.sampleRate;
    plot(x,y,'k-'); hold on
    
    set(h,'UserData',g);
else
    g = get(gcf,'UserData');
    s = g.s; T = g.T;
    for j = 1:length(varargin);
        argString = varargin{j};
        switch argString
            case 'next'
                if g.sweepNum < length(s)
                    g.sweepNum = g.sweepNum + 1;
                end
            case 'last'
                if g.sweepNum > 1
                    g.sweepNum = g.sweepNum - 1;
                end
            case 'jumpToSweep'
                if isempty(g.trialList)
                    nsweeps = s.length;
                    g.trialList = cell(1,nsweeps);
                    for k=1:nsweeps
                        g.trialList{k} = [int2str(k) ': trialNum=' int2str(s.trialNums(k))];
                    end
                end
                [selection,ok]=listdlg('PromptString','Select a sweep:','ListString',...
                    g.trialList,'SelectionMode','single');
                if ~isempty(selection) && ok==1
                    g.sweepNum = selection;
                end
            case 'rawSignal'
                g.displayType = 'rawSignal';
            case 'highPassFilteredSignal'
                g.displayType = 'highPassFilteredSignal';
            case 'LFP'
                g.displayType = 'LFP';
            case 'bitCode'
                g.displayType = 'bitCode';
            case 'none'
                g.displayTypeMinor = 'none';
            case 'spikeTimes'
                g.displayTypeMinor = 'spikeTimes';
            case 'spikeTimesPlusThreshold'
                g.displayTypeMinor = 'spikeTimesPlusThreshold';
            otherwise
                error('Invalid string argument.')
        end
    end
    
    subplot(2,1,1)
    cla; ms=8; lw=0.5;
    switch g.displayType
        case 'rawSignal'
            y=s.sweeps{g.sweepNum}.rawSignal;
            x=(1:length(y))/s.sweeps{g.sweepNum}.sampleRate;
            plot(x,y,'k-')
            xlabel('Sec');  ylabel('mV')
        case 'highPassFilteredSignal'
            y=s.sweeps{g.sweepNum}.highPassFilteredSignal;
            x=(1:length(y))/s.sweeps{g.sweepNum}.sampleRate;
            plot(x,y,'k-')
            xlabel('Sec');  ylabel('mV')
        case 'LFP'
            y=s.sweeps{g.sweepNum}.LFP;
            sampleRate = s.sweeps{g.sweepNum}.sampleRate;
            x=(1:length(y))/sampleRate;
            
           %---
%             bandPassCutOffsInHz = [5 25]; 
%             W1 = bandPassCutOffsInHz(1) / (sampleRate/2);
%             W2 = bandPassCutOffsInHz(2) / (sampleRate/2);
%             [b,a]=butter(2,[W1 W2]);
%             y = filtfilt(b, a, y);
            %----
            
            plot(x,y,'k-')
            xlabel('Sec');  ylabel('mV')
        case 'bitCode'
            y=s.sweeps{g.sweepNum}.bitCode;
            x=(1:length(y))/s.sweeps{g.sweepNum}.sampleRate;
            plot(x,y,'k-')
            xlabel('Sec');  ylabel('V')
    end
    
    switch g.displayTypeMinor
        case 'none'
            % do nothing
        case 'spikeTimes'
            hold on
            spikeTimeInds = s.sweeps{g.sweepNum}.spikeTimes;
            plot(x(spikeTimeInds),y(spikeTimeInds),'ro','MarkerSize',ms)
        case 'spikeTimesPlusThreshold'
            hold on
            spikeTimeInds = s.sweeps{g.sweepNum}.spikeTimes;
            plot(x(spikeTimeInds),y(spikeTimeInds),'ro','MarkerSize',ms)
            xmin = min(x); xmax = max(x);
            thresh = s.sweeps{g.sweepNum}.spikeThreshold;
            line([xmin xmax],[thresh thresh],'LineStyle','-','LineWidth',lw,'Color','g')
    end
    
    subplot(2,1,2); cla               
    if ~isempty(T.trials{g.sweepNum}.whiskerTrial)
        plotSymString = {'k-','r-','b-','g-','y-','m-','c-'};
%         tid = T.trials{g.sweepNum}.whiskerTrial.trajectoryIDs;
        tid=1;
        numTid = length(tid);
        for j=1:numTid
            s = plotSymString{mod(j,length(plotSymString))};
            [y,x] = T.trials{g.sweepNum}.whiskerTrial.get_position(tid(j));
            
            %---
%             sampleRate=500;
%             bandPassCutOffsInHz = [5 25]; 
%             W1 = bandPassCutOffsInHz(1) / (sampleRate/2);
%             W2 = bandPassCutOffsInHz(2) / (sampleRate/2);
%             [b,a]=butter(2,[W1 W2]);
%             y = filtfilt(b, a, y);
            %----
            
            plot(x,y,s);
        end
    else
        subplot(2,1,2); text(.1, .5, 'No whisker data this trial')
        
    end
    
    
    
end
% title([s.cellNum s.cellCode ', ' int2str(g.sweepNum) '/' int2str(s.length) ...
%     ', trialNum=' int2str(s.trialNums(g.sweepNum)) ...
%     ', xsgFileNum=' int2str(s.sweeps{g.sweepNum}.xsgFileNum) '\newline displayType=' g.displayType])

set(gcf,'UserData',g);
end







