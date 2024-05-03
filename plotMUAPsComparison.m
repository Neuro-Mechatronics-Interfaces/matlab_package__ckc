function fig = plotMUAPsComparison(MUAPs1,MUAPs2,options)
%PLOTMUAPSCOMPARISON Plot comparison of two template MUAP waveforms.
%
% Syntax:
%   fig = plotMUAPsComparison(MUAPs1, MUAPs2, 'Name', value, ...);
%
% arguments
%     MUAPs1
%     MUAPs2
%     options.SampleRate (1,1) {mustBePositive} = 2000;
%     options.BackColor (1,3) double {mustBeInRange(options.BackColor, 0, 1)} = [1 1 1];
%     options.FrontColor (1,3) double {mustBeInRange(options.FrontColor, 0, 1)} = [0 0 0];
%     options.FontSize (1,1) = 10;
% end

arguments
    MUAPs1
    MUAPs2
    options.SampleRate (1,1) {mustBePositive} = 2000;
    options.BackColor (1,3) double {mustBeInRange(options.BackColor, 0, 1)} = [1 1 1];
    options.FrontColor (1,3) double {mustBeInRange(options.FrontColor, 0, 1)} = [0 0 0];
    options.FontSize (1,1) = 10;
end
YLimMin=10000;
YLimMax=0;

hold on;

rY=size(MUAPs1,1);
cY=size(MUAPs1,2);
    
for r=1:rY
    for c=1:cY
        if ~isempty(MUAPs1{r,c})
            T=floor(size(MUAPs1{r,c},2)/2);
        end
    end
end
    
tmp = [];
for r=1:rY
    for c=1:cY
        mMUAPs =[];
        if isempty(MUAPs1{r,c})
            tmp(end+1) = 0;
        else
            for k1=1:size(MUAPs1{r,c},1)
                [mMUAPs(k1), ~]=max(abs([MUAPs1{r,c}(k1,:) MUAPs2{r,c}(k1,:)]));
            end
            tmp(end+1) = 1.0*max(mMUAPs)+eps;
        end
        tmp(end)=tmp(end)+eps;
    end
end
for r=1:rY
    dh(r) = 1.5*max(tmp);
end
yLineHeight = round(max(tmp)*100000)/100;
if yLineHeight > 10
    yLineHeight = floor(max(tmp)/10)*10;
end
xLineWidth = round(0.02*options.SampleRate);


tmpVec1 = []; tmpVec2 = [];
for r=1:size(MUAPs1,1)
    for c=1:size(MUAPs1,2)
        tmpVec1 = [tmpVec1 MUAPs1{r,c}];
        tmpVec2 = [tmpVec2 MUAPs2{r,c}];
    end
end
tcc = xcorr(tmpVec1,tmpVec2,20,'coeff');
[~,d] = max(abs(tcc));
d=d-20-1;
   
for r=1:rY
    for c=1:cY
        if ~isempty(MUAPs1{r,c})
            
            mMUAPs =[];
            if size(MUAPs1{r,c},1) > 1
                for k1=1:size(MUAPs1{r,c},1)
                    plot((c-1)*T*1.1+[1:T],MUAPs1{r,c}(k1,:).'+sum(dh(1:r)),'Color',[0.4,0.4,0.8]);
                    plot((c-1)*T*1.1+[1:T],MUAPs2{r,c}(k1,:).'+sum(dh(1:r)),'Color',[0.8,0.4,0.4]);
                end
            end
            xx = (c-1)*(2*T+1)*1.1+[1:(2*T+1)];
            if d>0
                plot(xx(1:end-d+1),mean(MUAPs1{r,c}(d:end),1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.1,0.1,0.8]);
                plot(xx(1:end-d+1),mean(MUAPs2{r,c}(1:end-d+1),1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.8,0.1,0.1]);
                %plot(xx(1:end-d+1),mean(MUAPs1{r,c}(d:end),1).'-mean(MUAPs2{r,c}(1:end-d+1),1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.1,0.8,0.1]);
            elseif d<0
                plot(xx(1:end+d+1),mean(MUAPs1{r,c}(1:end+d+1),1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.1,0.1,0.8]);       
                plot(xx(1:end+d+1),mean(MUAPs2{r,c}(-d:end),1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.8,0.1,0.1]);
                %plot(xx(1:end+d+1),mean(MUAPs1{r,c}(1:end+d+1),1).'-mean(MUAPs2{r,c}(-d:end),1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.1,0.8,0.1]);   
            else
                plot(xx,mean(MUAPs1{r,c},1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.1,0.1,0.8]);
                plot(xx,mean(MUAPs2{r,c},1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.8,0.1,0.1]);
                %plot(xx,mean(MUAPs1{r,c},1).'-mean(MUAPs2{r,c},1).'+sum(dh(1:r)),'LineWidth',0.5,'Color',[0.1,0.8,0.1]);
            end
            %plot3((c-1)*(2*T+1)*1.1+[1:(2*T+1)],sum(dh(1:r))+0*[1:(2*T+1)],mean(MUAPs1{r,c},1).','LineWidth',0.5,'Color',col{1}); view(-15,75);
            [mMUAPs(k1), ~]=max(abs([mean(MUAPs1{r,c},1) mean(MUAPs2{r,c},1)]));
        end
        
        YLimMax=max([YLimMax,max(mMUAPs)+sum(dh(1:r))]);
        YLimMin=min([YLimMin,min(mMUAPs)+sum(dh(1:r))]);
    end
end
    
plot((1)*(2*T+1)*1.1+[0,0],-dh(1)+[0,yLineHeight],'Color',options.FrontColor,'LineWidth',1);
plot((1)*(2*T+1)*1.1-[-T/3,+T/3],-dh(1) + [0,0],'Color',options.FrontColor,'LineWidth',1);
plot((1)*(2*T+1)*1.1-[-T/3,+T/3],-dh(1) + [yLineHeight,yLineHeight],'Color',options.FrontColor,'LineWidth',1);
text((1)*(2*T+1)*1.15+T/3,-dh(1)+yLineHeight/2,[num2str(yLineHeight) ' \muV'],'VerticalAlignment','middle','HorizontalAlignment','left','Color',options.FrontColor,'FontSize',options.FontSize);

axis tight;
hAx=gca;
set(hAx,'Color',options.BackColor);
set(hAx,'XColor',options.BackColor);
set(hAx,'YColor',options.BackColor);
set(hAx,'ZColor',options.BackColor);
set(hAx,'XTickLabel',[]);
set(hAx,'YTickLabel',[]);
title(d);
fig = gcf;
set(fig,'Color',options.BackColor);
end