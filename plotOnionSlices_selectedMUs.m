function [cum_IDR ,cum_EpochTime]=plotOnionSlices_selectedMUs(cMUPulses,fsamp,Force,BackColor,FrontColor)

MyColor=[ 
     0 0 1;...
     1 0 0;...
     0 1 0];

IDR_cutOff_freqency = 2.0; % cut-off frequency of low-pass fitler (in Hz)
fsamp2 = 100; % resampling frequency


[idr_B,idr_A] = butter (4, IDR_cutOff_freqency/fsamp2*2); % low-pass filter for instantaneous discharge rates [Hz]
[idr_B_short,idr_A_short] = butter (2, IDR_cutOff_freqency/fsamp2*2); % low-pass filter for instantaneous discharge rates [Hz]
outlierHTrsh=500/1000; %in s - high IPI threshold; for derecruitment of motor unit
outlierLTrsh=20/1000; % in s - low IPI threshold; for doublets

maxFir = -inf;
fir = {}; 

tmpForce = Force;
% [fB,fA] = butter (4, 2/fsamp*2); % low-pass filter
% for k1 = 1:size(Force,1)
%     if min(Force(k1,:))<0 & max(Force(k1,:)) < 0
%         Force(k1,:) = abs(Force(k1,:));
%     elseif min(Force(k1,:))<0
%         Force(k1,:) = Force(k1,:) - min(Force(k1,:));
%     end
%     tmpForce(k1,:) = filtfilt(fB,fA,Force(k1,:));
% end


for k1=1:length(cMUPulses); 
    if ~isempty(cMUPulses{k1})
        tmpMUPulses=sort(cMUPulses{k1})/fsamp;                        
        
        % instantaneous discharge rates and IPI variability plot
        tmpDichargPattern = diff(tmpMUPulses);
        tmpTime = tmpMUPulses(2:end);
        tmpInd=find(tmpDichargPattern > outlierLTrsh); %remove doublets
        tmpDichargPattern = tmpDichargPattern(tmpInd);
        tmpTime = tmpTime(tmpInd);
        
        tmpEpochs=find(tmpDichargPattern > outlierHTrsh); % find different MU recruitment epochs (i.e., epochs of MU activity)
        tmpEpochs=[0 tmpEpochs length(tmpDichargPattern)+1];
        
        for k2 = 1:length(tmpEpochs)-1
            tmpEpochTime = tmpTime(tmpEpochs(k2)+1:tmpEpochs(k2+1)-1);        
            %calculates mean IDR
            if length(tmpEpochTime) > 20; % filter epochs with more than 7 firings only
                tmpDichargPattern2 = interp1(tmpEpochTime,tmpDichargPattern(tmpEpochs(k2)+1:tmpEpochs(k2+1)-1),[tmpEpochTime(1):1/fsamp2:tmpEpochTime(end)]);
                tmpIPI = filtfilt(idr_B_short,idr_A_short,tmpDichargPattern2);
                %tmpIDR = 1./interp1([tmpEpochTime(1):1/fsamp2:tmpEpochTime(end)-1/fsamp2, tmpEpochTime(end)],tmpIPI,tmpEpochTime);
                tmpIDR = 1./tmpIPI;
                tmpEpochTime = [tmpEpochTime(1):1/fsamp2:tmpEpochTime(end)];
            elseif length(tmpEpochTime) > 10; % filter epochs with more than 7 firings only
                tmpDichargPattern2 = interp1(tmpEpochTime,tmpDichargPattern(tmpEpochs(k2)+1:tmpEpochs(k2+1)-1),[tmpEpochTime(1):1/fsamp2:tmpEpochTime(end)]);
                tmpIPI = filtfilt(idr_B_short,idr_A_short,tmpDichargPattern2);
                %tmpIDR = 1./interp1([tmpEpochTime(1):1/fsamp2:tmpEpochTime(end)-1/fsamp2, tmpEpochTime(end)],tmpIPI,tmpEpochTime); % the last value gets lost due to interpolation
                tmpIDR = 1./tmpIPI;
                tmpEpochTime = [tmpEpochTime(1):1/fsamp2:tmpEpochTime(end)];
            elseif length(tmpEpochTime) > 0 % for shorter epochs use the avarage IDR value
                tmpIDR = [];
                for k3 = tmpEpochs(k2)+1:tmpEpochs(k2)+length(tmpEpochTime)-3;
                    tmpIDR(end+1) = 1/mean(tmpDichargPattern(k3:k3+3));
                end
                tmpIDR(end+1:length(tmpEpochTime)) = 1/mean(tmpDichargPattern(max([tmpEpochs(k2)+1,tmpEpochs(k2+1)-4]):tmpEpochs(k2+1)-1));
                %tmpIDR = 1/mean(tmpDichargPattern(tmpEpochs(k2)+1:tmpEpochs(k2+1)-1))*ones(size(tmpEpochTime));
            else
                continue;
            end
                       
            cum_IDR{k1}{k2} = tmpIDR;            
            cum_EpochTime{k1}{k2} = tmpEpochTime;  
            maxFir = max([maxFir tmpIDR]);
        end

    end
end

hold on;
forceStr{1} = 'force x';
forceStr{2} = 'force y';
for k1 = 1:size(Force,1)
    plot([1:length(tmpForce(k1,:))]/fsamp,tmpForce(k1,:)/max([max(abs(tmpForce(:))),10^-10])*maxFir,'LineWidth',4,'Color',[0.7 0.7 0.7]+(1-k1)*0.3); 
    [m,i] = max(tmpForce(k1,:));
    line([i/fsamp,i/fsamp+1],[maxFir,maxFir+1],'LineWidth',2,'Color',FrontColor);
    ht = text(i/fsamp+1.1,maxFir+1,forceStr{k1});
    set(ht,'Color',FrontColor);
    set(ht,'FontSize',18);
    set(ht,'HorizontalAlignment','left');
end

for k1=1:length(cum_IDR);
    for k2 = 1:length(cum_IDR{k1})
        plot(cum_EpochTime{k1}{k2},cum_IDR{k1}{k2},'LineWidth',1,'Color',MyColor(mod(k1-1,size(MyColor,1))+1,:));
    end

    if length(cum_EpochTime{k1})>0
        ht = text(cum_EpochTime{k1}{k2}(end)+0.1,cum_IDR{k1}{k2}(end),['MU ' num2str(k1)]);
        set(ht,'Color',FrontColor);
        set(ht,'FontSize',10);
    end
end

xlabel('time (s)','FontSize',10);
ylabel('smoothed discharge rate (pps)','FontSize',10);

hAx=gca;   
set(hAx,'Color',BackColor);
set(hAx,'XColor',FrontColor);
set(hAx,'YColor',FrontColor);
set(hAx,'ZColor',FrontColor);   
set(hAx,'FontSize',10);   