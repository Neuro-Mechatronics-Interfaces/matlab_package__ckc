function [startRef, stopRef, startSample, stopSample] = findMostStableRegion(ref_signal,IntervalLengthInSamples,fsamp)

startRef = 1;
D = max(abs(ref_signal)) - min(abs(ref_signal));
stopRef = find(ref_signal > ref_signal(1)+D/10);
stopRef = min(stopRef);
stopRef = max(fsamp/2,stopRef-fsamp/2);

stopSig = find(ref_signal > ref_signal(1)+D/10);
stopSig = max(stopSig);

varRef = zeros(1,length(ref_signal) - IntervalLengthInSamples);
for startSample = stopRef:stopSig - IntervalLengthInSamples
    varRef(startSample) = var(ref_signal(startSample:startSample+IntervalLengthInSamples));  % var
    % varRef(startSample) = std(ref_signal(startSample:startSample+IntervalLengthInSamples))/mean(ref_signal(startSample:startSample+IntervalLengthInSamples));  %Cov
end
varRef(varRef==0) = max(varRef);

[~,startSample] = min(varRef);
stopSample = startSample + IntervalLengthInSamples;

end