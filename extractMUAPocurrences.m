function MUAP=extractMUAPocurrences(Y,MU_IND,WinSemiLength,winMode,TransSec)
%EXTRACTMUAPOCURRENCES  Routine for spike triggered extraction of MUAPs from multichannel measurements
% MUAP=extractMUAPocurrences(Y,MU_IND,WinSemiLength,winMode,TransSec)
%
% INPUTS:
%       Y - row-array of all the measurements (each row corresponds to one measurement)
%       MU_IND - MU's discharge times (in samples)
%       WinSemiLength - semilength of window averaging function
%       winMode - window type (optional): if not specified, rectangular window is used
%       TransSec - optional parameter for the Tukey window
% OUTPUTS:
%       MUAP - cell structure containing the row-array of all MUAP ocurrences for the specific measurements
%              MUAP{i} for i-th measurement (i.e. i-th row in Y).

if iscell(Y)
    sigLength = 0;
    for r = 1:size(Y,1);
        for c = 1:size(Y,2);
            if ~isempty(Y{r,c})
                sigLength = max([sigLength,length(Y{r,c})]);
            end
        end
    end
    MU_IND=MU_IND(find( MU_IND>WinSemiLength+1 & MU_IND<sigLength-WinSemiLength-1));
else
    MU_IND=MU_IND(find( MU_IND>WinSemiLength+1 & MU_IND<length(Y)-WinSemiLength-1));
end
rY=size(Y,1);

if (nargin > 3)
    WinSemiLength=WinSemiLength+mod(WinSemiLength,2); % make WinSemiLength even
    if (strcmp('BLACKMAN',upper(winMode)))
        wind=blackman(WinSemiLength+1)'; % Use Blackman window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('BLACKMAN-HARRIS',upper(winMode)))
        wind=blackmanharris(WinSemiLength+1)'; % Use Blackman-Harris window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('NUTTALLWIN',upper(winMode)))
        wind=nuttallwin(WinSemiLength+1)'; % Use nutall window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('HANN',upper(winMode)))
        wind=hann(WinSemiLength+1)'; % Use Hanning window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('CHEBWIN',upper(winMode)))
        wind=chebwin(WinSemiLength+1,100)'; % Use Chebyshev window
        wind=[wind(1:WinSemiLength/2) ones(1,WinSemiLength) wind(WinSemiLength/2+1:end)];
    elseif (strcmp('TUKEYWIN',upper(winMode)))
        wind=tukeywin(2*WinSemiLength+1,TransSec)'; % Use Tukey window
    else
        disp('ERROR: window not found');
        return;
    end
else
    wind=ones(1,2*WinSemiLength+1);
end

if isempty(MU_IND)
    for k=1:rY
        MUAP{k}=zeros(1,2*WinSemiLength+1);
    end
else
    
    if ~iscell(Y)
        for k=1:rY
            tmp = ckc.cutMUAP(MU_IND,WinSemiLength,Y(k,:));
            for k2 = 1:size(tmp,1)
                MUAP{k}(k2,:) = wind.*tmp(k2,:);
            end
        end
    else
        for r=1:size(Y,1)
            for c=1:size(Y,2)
                if ~isempty(Y{r,c})
                    tmp = ckc.cutMUAP(MU_IND,WinSemiLength,Y{r,c});
                    for k2 = 1:size(tmp,1)
                        MUAP{r,c}(k2,:) = wind.*tmp(k2,:);
                    end
                end
            end
        end
    end
    
end