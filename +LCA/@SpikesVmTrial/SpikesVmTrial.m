classdef SpikesVmTrial < LCA.SpikesTrial
    %
    %
    %
    %
   properties

        rawSignal=[];
        AOM=[];
        FP=[];
    end



    methods (Access = public)
        function obj = SpikesVmTrial(trial_num, cell_num, cell_code, xsg_file_num,...
                sample_rate, sweep_length_in_samples, spike_times, rawSignal, AOM, FP)
            %
            %  function obj = SpikesTrial(trial_num, cell_num, cell_code, xsg_file_num,...
            %    sample_rate, spike_times)
            %
            %
            
             obj = obj@LCA.SpikesTrial(trial_num, cell_num, cell_code, xsg_file_num,...
                sample_rate, sweep_length_in_samples, spike_times);
             obj.rawSignal = rawSignal;
             obj.AOM=AOM;
             if nargin<10
                 obj.FP=[];
             else
                 obj.FP=FP;
             end;
        end

   
    end
end

