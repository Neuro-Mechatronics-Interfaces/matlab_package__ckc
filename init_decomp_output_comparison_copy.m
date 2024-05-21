function target_hat = init_decomp_output_comparison_copy(target, IPTs, options)
%INIT_DECOMP_OUTPUT_COMPARISON_COPY  Initialize a copy of the target output for computing synchronization index etc.
arguments
    target
    IPTs
    options.MinPeakHeight (1,1) double = 0.35;
end
target_hat = target;
n = size(IPTs,1);
target_hat.IPTs = IPTs;
target_hat.MUPulses = cell(1,n);
for ii = 1:n
    [~,target_hat.MUPulses{ii}] = findpeaks(IPTs(ii,:),'MinPeakHeight',options.MinPeakHeight);
end
end