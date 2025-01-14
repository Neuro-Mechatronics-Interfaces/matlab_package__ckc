function fig = plotMUAPs(MUAP,fsamp, options)
%PLOTMUAPS  Plot MUAP template waveforms
%
% Syntax:
%   fig = ckc.plotMUAPs(MUAP, fsamp, 'Name', value, ...);
%
% Example:
%   data = load(fullfile(F(iF).folder,F(iF).name));
%   label = strrep(data.description,": ","--");
%   Description{iF} = repmat(string(data.description), numel(data.MUPulses), 1);
%   [PNR{iF},MUID_Label] = batch_compute_PNR_and_Label(data);
%   MyColor=COLORMAP_BASE(end : -floor(size(COLORMAP_BASE,1)/max([length(data.MUPulses),1])) :1,:);
%   for iMUAP = 1:numel(data.MUPulses)
%       [MUAP,MUAPrms,MUAPp2p,MUAPduration,AllMUAPs] = ckc.extractMUAPs(data.SIG,data.MUPulses{iMUAP},data.fsamp);
%       fig = ckc.plotMUAPs(AllMUAPs, data.fsamp, ...
%           'FrontColor', MyColor(iMUAP,:), ...
%           'Title', MUID_Label{iMUAP}, ...
%           'Subtitle', sprintf('(N = %d Firings)', size(AllMUAPs{1,1},1)));
%       utils.save_figure(fig, fullfile(OUTPUT_ROOT,"Templates",label), sprintf("MUAP-%02d-Template", iMUAP), 'ExportAs', {'.png'});
%   end
arguments
    MUAP
    fsamp (1,1) {mustBePositive}
    options.BackColor  = [1, 1, 1];
    options.FrontColor = [0, 0, 0];
    options.IndividualTraceColor = [0.65 0.65 0.65];
    options.OverlayIndividual = true;
    options.MaxIndividualTraces = 25;
    options.Title {mustBeTextScalar} = "";
    options.Subtitle {mustBeTextScalar} = "";
    options.FontName = 'Tahoma';
    options.FontSize = 14;
    options.YScale = 0.75;
    options.YBarHeight = 100; % uV
end

YLimMin=10000;
YLimMax=0;

fig=figure('Name', 'MUAP Template Waveforms', ...
    'WindowState', 'maximized', ...
    'PaperPositionMode', 'auto', ...
    'Color', options.BackColor);

rY=size(MUAP,1);
cY=size(MUAP,2);
figInRow = cY;


for r=1:rY
    for c=1:cY
        if ~isempty(MUAP{r,c})
            T=floor(size(MUAP{r,c},2)/2);
        end
    end
end

tmp = [];
for r=1:rY
    for c=1:cY
        mMUAP1 =[];
        if isempty(MUAP{r,c})
            tmp(end+1) = 0;
        else
            for k1=1:size(MUAP{r,c},1)
                [mMUAP1(k1), iMUAP1]=max(abs(MUAP{r,c}(k1,:)));
            end
            tmp(end+1) = 1.0*max(mMUAP1)+eps;
        end
        tmp(end)=tmp(end)+eps;
    end
end

dh = options.YScale*max(tmp);

if isempty(options.YBarHeight)
    yLineHeight = round(dh*1000)/1000;
    if yLineHeight > 10
        yLineHeight = floor(max(tmp)/10)*10;
    end
else
    yLineHeight = options.YBarHeight;
end

xLineWidth = 10/(2*T/fsamp*1000)/1.05;

for k=1:figInRow
    All_iMUAP1{k} = [];
end

hold on;
vec = reshape(randsample(size(MUAP{1,1},1),min(size(MUAP{1,1},1),options.MaxIndividualTraces),false),1,[]);
for r=1:rY
    for c=1:cY
        if ~isempty(MUAP{r,c})
            %subplot(1,figInRow,c); hold on;
            mMUAP1 =[];
            
            if options.OverlayIndividual
                for k1=vec
                    plot(c+(1:size(MUAP{r,c},2))/(size(MUAP{r,c},2)*1.05),MUAP{r,c}(k1,:).'+dh*r,'Color',options.IndividualTraceColor);
                end
            end
            plot(c+(1:size(MUAP{r,c},2))/(size(MUAP{r,c},2)*1.05), mean(MUAP{r,c},1).'+dh*r,'LineWidth',1.5,'Color',options.FrontColor); %col{2});
        end

        YLimMax=max([YLimMax,max(mMUAP1)+dh*r]);
        YLimMin=min([YLimMin,min(mMUAP1)+dh*r]);
    end
end
    
hAx=get(fig,'CurrentAxes'); 
set(hAx,'YTick',[]);
set(hAx,'XTick',[]);
set(hAx,'Units','normalized');
set(hAx,'Position',[0.02,0.02,0.96,0.9]);

plot(cY+xLineWidth*[1,2],-dh/2+ [0,0],'Color',options.FrontColor,'LineWidth',2);
plot(cY+xLineWidth*[1,1],-dh/2+ [-yLineHeight/10,+yLineHeight/10],'Color',options.FrontColor,'LineWidth',2);
plot(cY+xLineWidth*[2,2],-dh/2+ [-yLineHeight/10,+yLineHeight/10],'Color',options.FrontColor,'LineWidth',2);
text(cY+xLineWidth*1.5,-dh/2-yLineHeight/10-10,'10 ms','VerticalAlignment','middle','HorizontalAlignment','center','Color',options.FrontColor,'FontName','Tahoma','FontSize',options.FontSize);

plot(cY+xLineWidth*[3,3],-dh/2+[0,yLineHeight],'Color',options.FrontColor,'LineWidth',2);
plot(cY+xLineWidth*3+[-xLineWidth/3,+xLineWidth/3],-dh/2 + [0,0],'Color',options.FrontColor,'LineWidth',2);
plot(cY+xLineWidth*3+[-xLineWidth/3,+xLineWidth/3],-dh/2 + [yLineHeight,yLineHeight],'Color',options.FrontColor,'LineWidth',2);
text(cY+xLineWidth*3+1.5*xLineWidth/3,-dh/2+yLineHeight/2,[num2str(yLineHeight) ' \muV'],'VerticalAlignment','middle','HorizontalAlignment','left','Color',options.FrontColor,'FontName','Tahoma','FontSize',options.FontSize);

set(hAx,'Color',options.BackColor);
set(hAx,'XColor',options.BackColor);
set(hAx,'YColor',options.BackColor);
set(hAx,'ZColor',options.BackColor);

if strlength(options.Title) > 0
    title(hAx, options.Title, 'FontName', options.FontName, 'Color', options.FrontColor);
end

if strlength(options.Subtitle) > 0
    subtitle(hAx, options.Subtitle, 'FontName', options.FontName, 'Color', [0.65 0.65 0.65]);
end

end