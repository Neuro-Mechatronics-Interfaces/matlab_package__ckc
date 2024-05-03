function PNR = calcPNR(MUPulses,IPT,fsamp)
%COMPUTE_PNR  Calculates the pulse-to-noise ratio (PNR) 
%
% Syntax:
%   PNR = ckc.calcPNR(MUPulses, IPT, fsamp);

arguments
    MUPulses (1,:) double
    IPT (1,:) double
    fsamp (1,1) {mustBePositive}
end

if ~isempty(MUPulses)
    
    surrogatePulseSamples = [];
    for t = -round(3/2048*fsamp):round(3/2048*fsamp)
        surrogatePulseSamples = [surrogatePulseSamples MUPulses+t]; %#ok<AGROW>
    end
    
    nInd = setdiff(MUPulses(1):MUPulses(end),surrogatePulseSamples);
    IPT = IPT/mean(IPT(MUPulses));
    
    tmpn = IPT(nInd);
    tmpn = tmpn(~isnan(tmpn));
    tmpn = tmpn(tmpn>=0);
    PNR = round(10*10*log10(  mean(IPT(MUPulses).^2) / mean(tmpn.^2)))/10;
else
    PNR = 0;
end

end