function fig = plotIDR(MUPulses,MU_ID,ref_signal,fsamp,options)
%PLOT_IDR Plots MU instantaneous discharge rate.
% 
% Syntax: 
%   fig = plotIDR(MUPulses,MU_ID,ref_signal,'Name',value,...)
%
% Input:
%    MUPulses - cell structure containing the discharge patterns of all MUs, e.g. MUPulses{i} contains firing instants of the i-th MU (in samples).
%    MU_ID - IDs of motor units to display
%    ref_signal - array containing the measured contraction ref_signal
%    fsamp - sampling frequency

arguments
    MUPulses (1,:) cell
    MU_ID    (1,:) cell
    ref_signal (:,:) double
    fsamp (1,1) {mustBePositive}
    options.BackColor (1,3) double {mustBeInRange(options.BackColor, 0, 1)} = [1 1 1];
    options.FrontColor (1,3) double {mustBeInRange(options.FrontColor, 0, 1)} = [0 0 0];
    options.ColorMap (256,3) double {mustBeInRange(options.ColorMap,0,1)} = jet(256);
    options.Title {mustBeTextScalar} = "";
    options.Subtitle {mustBeTextScalar} = "";
    options.YLabelLeft {mustBeTextScalar} = "Motor Unit";
    options.YLabelRight {mustBeTextScalar} = "Instantaneous Discharge Rate (pps)";
    options.FontSize (1,1) {mustBePositive} = 12;
    options.FontName {mustBeTextScalar} = "Tahoma";
    options.XLabel {mustBeTextScalar} = "Time (s)";
    options.RefSignalWidth (1,1) {mustBePositive} = 3;
    options.RefSignalBaseColor (1,3) {mustBeInRange(options.RefSignalBaseColor,0,1)} = [0.7 0.7 0.7];
end

fig=figure('Name','MUAP IDR Overlay');

MyColor=colormap(options.ColorMap);
MyColor=MyColor(end : -floor(size(MyColor,1)/max([length(MUPulses),1])) :1,:);

outlierHTrsh=500*fsamp/1000; %500; % outlier threshold in ms
outlierLTrsh=15*fsamp/1000; % 10; % outlier threshold in ms
d=30; % options.FontSize0;
dH=20; % 80;

LLim=1;
HLim=length(ref_signal);

tmpref_signal = ref_signal;

% Logic to remove scatter points outside the time-range of the reference
% signal. 
r=length(MUPulses);
tmp=cell(1,r);
tmp2=cell(1,r);
lenTmp=nan(1,r);
for k=1:r
    tmp2{k}=MUPulses{k};
    MU=zeros(1,HLim); MU(tmp2{k})=1;
    tmp2{k}=find(MU(LLim:HLim));
    tmp{k}=(tmp2{k}(2:end)-tmp2{k}(1:end-1))/fsamp*1000;
    tmpInd=find(outlierLTrsh<tmp{k} & tmp{k}<outlierHTrsh);
    tmp{k}=1000./tmp{k}(tmpInd);
    tmp2{k}=tmp2{k}(tmpInd);
    lenTmp(k)=length(tmp{k});
end

% Create the "left" and "right" axes
newAxes = axes(fig,'NextPlot','add','Position',[0.06 80/max(min(r*150,970),600) 0.80 0.9-80/max(min(r*150,970),600)]); %,'Visible','off');
hAx =     axes(fig,'NextPlot','add','Position',[0.06 80/max(min(r*150,970),600) 0.80 0.9-80/max(min(r*150,970),600)]);

% Overlay reference signals first
for k1 = 1:size(ref_signal,1)
    tmpref_signal(k1,:) = (r+0.5)*(d+dH)*abs(tmpref_signal(k1,LLim:HLim))/max([max(abs(tmpref_signal(:))),10^-10]);
    plot(hAx,tmpref_signal(k1,:), ...
        'LineWidth',options.RefSignalWidth,'Color',options.RefSignalBaseColor+(1-k1)*0.3);
end

% Iterate over MUAPs, adding the instantaneous firings as individual points
for k=1:r
    if length(tmp{k})>1
        plot(tmp2{k}(1:end),(k-1)*(d+dH)+tmp{k},'o', ...
            'MarkerSize',4, ...
            'LineWidth',0.1, ...
            'Color',MyColor(mod(k-1,size(MyColor,1))+1,:), ...
            'MarkerFaceColor',MyColor(mod(k-1,size(MyColor,1))+1,:));
    end
end

YTickValue = sort([d/3:(dH+d):(r)*(dH+d) 2*d/3:(dH+d):r*(dH+d) d:(dH+d):r*(dH+d)]);
YTickLabel=cell(1,3);
YTickLabel{1}=sprintf('%d', round(YTickValue(1))); 
YTickLabel{2}='';
YTickLabel{3}=sprintf('%d', round(YTickValue(3)));

set(hAx,'YAxisLocation','right', ...
        'YLim',[0,(r+0.5)*(d+dH)+dH], ...
        'YTick',YTickValue, ...
        'YTickLabel',YTickLabel, ...
        'XLim',[1,HLim-LLim+1]);

if(HLim/fsamp > 20)
    step=5;
else
    step=1;
end
XTickValues = (1:step*fsamp:HLim-LLim+1)+1;
XTickLabels=round((LLim+XTickValues)/fsamp*10)/10;

set(hAx,'XTick',XTickValues, ...
        'XTickLabel',XTickLabels, ...
        'FontSize',options.FontSize);
xlabel(hAx,options.XLabel, 'FontSize', options.FontSize, 'FontName', options.FontName);
ylabel(hAx,options.YLabelRight, 'FontSize', options.FontSize, 'FontName', options.FontName);


% Label the MUAPs axes using the corresponding colors
set(fig,'CurrentAxes',newAxes)
MUAP_Label_Tick_Values = sort((dH+d)/2:(dH+d):(r-0.5)*(dH+d));
MU_Labels = cell(size(MU_ID));
for ii = 1:numel(MU_ID)
    MU_Labels{ii} = sprintf('\\bf\\color[rgb]{%04.2f,%04.2f,%04.2f}%s',...,
        round(MyColor(ii,:),2),MU_ID{ii});
end

yline(newAxes,MUAP_Label_Tick_Values+(d+dH)/2, ...
    'LineStyle',':','Color',[0.65 0.65 0.65],'LineWidth',1);
set(newAxes,'XTick',[], ...
            'xTickLabel',[], ...
            'YLim',[0,(r+0.5)*(d+dH)+dH], ...
            'YTick',MUAP_Label_Tick_Values, ...
            'YTickLabel',MU_Labels, ...
            'FontSize',options.FontSize);
ylabel(newAxes, options.YLabelLeft, 'FontSize', options.FontSize);

set(fig,'Color',options.BackColor, ...
      'WindowState','maximized', ...
      'PaperPositionMode','auto');
set(hAx,'Color','None', ...
    'XColor',options.FrontColor, ...
    'YColor',options.FrontColor, ...
    'ZColor',options.FrontColor, ...
    'FontName','Tahoma');

set(newAxes,'Color','None', ...
    'XColor',options.FrontColor, ...
    'YColor',options.FrontColor, ...
    'ZColor',options.FrontColor, ...
    'FontName','Tahoma');

if strlength(options.Title) > 0
    title(newAxes,options.Title,'FontName','Tahoma','Color',options.FrontColor);
end

if strlength(options.Subtitle) > 0
    subtitle(newAxes,options.Subtitle,'FontName','Tahoma','Color',[0.65 0.65 0.65]);
end
end