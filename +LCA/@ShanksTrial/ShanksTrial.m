classdef ShanksTrial < handle
    %
    %
    %
    %
    % DHO, 5/08.
    %
    %
    %
    %
    properties
        shanks      = []; %: 108
        clustData   = [];
        multis      = [];
        useFlag     = 1;
    end

    properties (Dependent=true)
        trialNum
        numClusts
        sampleRate
        sweepLengthInSamples
        spikeRateInHz
    end

    methods (Access = public)
        function obj = ShanksTrial(shanks, clusts)

            if nargin > 0
                obj.shanks  = shanks;
                obj.clustData = clusts;
            end
        end

        
        function r = getSpikeTimes(obj, shankNum, cellNum)
            %
            % r is in seconds.
            %
            if  isempty(obj.clustData)
                r = [];

            else
                s = find(cellfun(@(x)find(x.cellNum == cellNum & x.shankNum == shankNum),obj));
                r = s ./ obj.sampleRate;
            end
        end
        

     
       
    end

    methods % Dependent property methods; cannot have attributes.
     
        function value = get.trialNum(obj)
            if ~isempty(obj.clustData)
                value = (obj.clustData{1}.trialNum);
            else
                value = [];
            end
        end
        
         function value = get.sampleRate(obj)
            if ~isempty(obj.clustData)
                value = (obj.clustData{1}.sampleRate);
            else
                value = [];
            end
         end
        
          function value = get.sweepLengthInSamples(obj)
            if ~isempty(obj.clustData)
                value = (obj.clustData{1}.sweepLengthInSamples);
            else
                value = [];
            end
        end
        
        function value = get.numClusts(obj)
            if  isempty(obj.clustData)
                value = [];
        
            else
                value = length(obj.clustData);
            end
        end
        
               
        
         function value = get.spikeRateInHz(obj)
            if  isempty(obj.clustData)
                value = [];
            else
                value = cellfun(@(x)x.spikeRateInHz,obj.clustData,'UniformOutput',0);

            end
        end

    end
end

