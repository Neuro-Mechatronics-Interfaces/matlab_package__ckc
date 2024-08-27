function [fig, MUPulses, IDR_Value] = plotIDR(MUPulses,MU_ID,ref_signal,fsamp,options)
%PLOT_IDR Plots MU instantaneous discharge rate.
% 
% Syntax: 
%   fig = ckc.plotIDR(MUPulses,MU_ID,ref_signal,'Name',value,...)
%   [fig, MUPulses, IDR_Value] = ckc.plotIDR(...);
%
% Input:
%    MUPulses - cell structure containing the discharge patterns of all MUs, e.g. MUPulses{i} contains firing instants of the i-th MU (in samples).
%    MU_ID - IDs of motor units to display
%    ref_signal - array containing the measured contraction ref_signal
%                   -> nRefChannels x nRefSamples
%    fsamp - sampling frequency
%
% Options:
%   BackColor (1,3) double {mustBeInRange(options.BackColor, 0, 1)} = [1 1 1];
%   FrontColor (1,3) double {mustBeInRange(options.FrontColor, 0, 1)} = [0 0 0];
%   ColorMap (256,3) double {mustBeInRange(options.ColorMap,0,1)} = jet(256);
%   Title {mustBeTextScalar} = "";
%   Subtitle {mustBeTextScalar} = "";
%   SortMethod {mustBeMember(options.SortMethod, {'None', 'Recruitment'})} = 'Recruitment';
%   YLabelLeft {mustBeTextScalar} = "Motor Unit";
%   YLabelRight {mustBeTextScalar} = "Instantaneous Discharge Rate (pps)";
%   FontSize (1,1) {mustBePositive} = 12;
%   FontName {mustBeTextScalar} = "Tahoma";
%   XLabel {mustBeTextScalar} = "Time (s)";
%   RateOutlierThreshold (1,2) = [3, 50]; % Pulses/sec lower and upper outlier limits
%   TimeOutlierThreshold (1,2) = [5, inf]; % Exclude samples before/after this number of seconds
%   HeightSpacing (1,1) double = 30; % "d"
%   HeightPadding (1,1) double = 20; % "dH"
%   RefSignalWidth (1,1) {mustBePositive} = 3;
%   RefSignalBaseColor (1,3) {mustBeInRange(options.RefSignalBaseColor,0,1)} = [0.7 0.7 0.7];
%
% Output:
%   fig - Figure handle of resulting figure
%   MUPulses - 1 x nMUAPs cell array of pulse instants for plotting IDRs
%                   -> Removes Outlier Exclusions based on IDR_Value; see
%                       'RateOutlierThreshold' option for details. Pulses
%                       can also be removed due to 'TimeOutlierThreshold'
%                       setting.
%   IDR_Value  - 1 x nMUAPs cell array of IDR values used in plot
%
% See also: Contents

arguments
    MUPulses (1,:) cell
    MU_ID    (1,:) cell
    ref_signal (:,:) double % nRefChannels x nRefSamples signal
    fsamp (1,1) {mustBePositive}
    options.BackColor (1,3) double {mustBeInRange(options.BackColor, 0, 1)} = [1 1 1];
    options.FrontColor (1,3) double {mustBeInRange(options.FrontColor, 0, 1)} = [0 0 0];
    options.ColorMap (256,3) double {mustBeInRange(options.ColorMap,0,1)} = jet(256);
    options.Title {mustBeTextScalar} = "";
    options.Subtitle {mustBeTextScalar} = "";
    options.ManualSortOrder = [];
    options.SortMethod {mustBeMember(options.SortMethod, {'None', 'Manual', 'Recruitment'})} = 'Recruitment';
    options.YLabelLeft {mustBeTextScalar} = "Motor Unit";
    options.YLabelRight {mustBeTextScalar} = "Instantaneous Discharge Rate (pps)";
    options.FontSize (1,1) {mustBePositive} = 12;
    options.FontName {mustBeTextScalar} = "Tahoma";
    options.XLabel {mustBeTextScalar} = "Time (s)";
    options.MinPulsesRequired (1,1) double {mustBeInteger, mustBePositive} = 20;
    options.RateOutlierThreshold (1,2) = [1, 100];
    options.TimeOutlierThreshold (1,2) = [1, inf]; % Exclude samples before/after this number of seconds
    options.HeightSpacing (1,1) double = 45; % "d"
    options.HeightPadding (1,1) double = 5; % "dH"
    options.RefSignalWidth (1,1) {mustBePositive} = 1.25;
    options.RefSignalBaseColor (1,3) {mustBeInRange(options.RefSignalBaseColor,0,1)} = [0.7 0.7 0.7];
end

fig=figure('Name','MUAP IDR Overlay');

MyColor=colormap(options.ColorMap);
MyColor=MyColor(end : -floor(size(MyColor,1)/max([length(MUPulses),1])) :1,:);

sampleIndexThreshold = options.TimeOutlierThreshold .* fsamp;
nSamplesRef = size(ref_signal,2);

% Calculate the actual IDR values
r=numel(MUPulses);
IDR_Value=cell(1,r);
for k=1:r
    valid_times  = (MUPulses{k} > sampleIndexThreshold(1)) & ...
                   (MUPulses{k} < sampleIndexThreshold(2)); 
    MUPulses{k} = MUPulses{k}(valid_times);
    IDR_Value{k} = fsamp./diff(MUPulses{k});
    valid_rates  = (IDR_Value{k} > options.RateOutlierThreshold(1)) & ...
                   (IDR_Value{k} < options.RateOutlierThreshold(2));
    IDR_Value{k} = IDR_Value{k}(valid_rates);
    MUPulses{k}  = MUPulses{k}(valid_rates);
end

% Create the "left" and "right" axes
MUAPLabelAxes_Left = axes(fig,'NextPlot','add','Position',[0.06 80/max(min(r*150,970),600) 0.80 0.9-80/max(min(r*150,970),600)]); %,'Visible','off');
IDRAxes_Right      = axes(fig,'NextPlot','add','Position',[0.06 80/max(min(r*150,970),600) 0.80 0.9-80/max(min(r*150,970),600)]);

% Overlay reference signals first
for k1 = 1:size(ref_signal,1)
    ref_signal(k1,:) = (r+0.5)*(options.HeightSpacing+options.HeightPadding)*abs(ref_signal(k1,:))/max([max(abs(ref_signal(:))),10^-10]);
    plot(IDRAxes_Right,ref_signal(k1,:), ...
        'LineWidth',options.RefSignalWidth,'Color',options.RefSignalBaseColor+(1-k1)*0.3);
end

switch options.SortMethod
    case 'None'
        sortVector = 1:r;
    case 'Manual'
        if numel(options.ManualSortOrder) == r
            sortVector = reshape(options.ManualSortOrder,1,r);
        else
            warning('ManualSortOrder must have exactly %d elements (corresponding to indexed elements of MUPulses input). Using "None" for SortMethod instead.', r);
            sortVector = 1:r;
        end
    case 'Recruitment'
        sortVector = ckc.sortByRecruitmentOrder(MUPulses);
    otherwise
        error("SortMethod '%s' not yet implemented!", options.SortMethod);
end

% Iterate over MUAPs, adding the instantaneous firings as individual points
for ii=1:r
    k = sortVector(ii);
    nPulses = numel(IDR_Value{k});
    if nPulses > options.MinPulsesRequired
        IDR_Offset = (ii-1)*(options.HeightSpacing+options.HeightPadding);
        line(IDRAxes_Right, ...
            MUPulses{k}, ...
            IDR_Value{k}+IDR_Offset, ...
            'Marker', 'o', ...
            'MarkerSize',4, ...
            'LineStyle', 'none', ...
            'MarkerEdgeColor',MyColor(ii,:), ...
            'MarkerFaceColor',MyColor(ii,:), ...
            'DisplayName', MU_ID{k});
    end
end

YTickValue = sort([options.HeightSpacing/3:(options.HeightPadding+options.HeightSpacing):(r)*(options.HeightPadding+options.HeightSpacing) 2*options.HeightSpacing/3:(options.HeightPadding+options.HeightSpacing):r*(options.HeightPadding+options.HeightSpacing) options.HeightSpacing:(options.HeightPadding+options.HeightSpacing):r*(options.HeightPadding+options.HeightSpacing)]);
YTickLabel=cell(1,3);
YTickLabel{1}=sprintf('%d', round(YTickValue(1))); 
YTickLabel{2}='';
YTickLabel{3}=sprintf('%d', round(YTickValue(3)));

set(IDRAxes_Right,'YAxisLocation','right', ...
        'YLim',[0,(r+0.5)*(options.HeightSpacing+options.HeightPadding)+options.HeightPadding], ...
        'YTick',YTickValue, ...
        'YTickLabel',YTickLabel, ...
        'XLim',[1,nSamplesRef]);

if(nSamplesRef/fsamp > 20)
    step=5;
else
    step=1;
end
XTickValues = (1:step*fsamp:nSamplesRef)+1;
XTickLabels=round(XTickValues/fsamp*10)/10;

set(IDRAxes_Right,'XTick',XTickValues, ...
        'XTickLabel',XTickLabels, ...
        'FontSize',options.FontSize);
xlabel(IDRAxes_Right,options.XLabel, 'FontSize', options.FontSize, 'FontName', options.FontName);
ylabel(IDRAxes_Right,options.YLabelRight, 'FontSize', options.FontSize, 'FontName', options.FontName);


% Label the MUAPs axes using the corresponding colors
set(fig,'CurrentAxes',MUAPLabelAxes_Left)
MUAP_Label_Tick_Values = sort((options.HeightPadding+options.HeightSpacing)/2:(options.HeightPadding+options.HeightSpacing):(r-0.5)*(options.HeightPadding+options.HeightSpacing));
MU_Labels = cell(size(MU_ID));
for ii = 1:numel(MU_ID)
    MU_Labels{ii} = sprintf('\\bf\\color[rgb]{%04.2f,%04.2f,%04.2f}%s',...,
        round(MyColor(ii,:),2),MU_ID{sortVector(ii)});
end

yline(MUAPLabelAxes_Left,MUAP_Label_Tick_Values+(options.HeightSpacing+options.HeightPadding)/2, ...
    'LineStyle',':','Color',[0.65 0.65 0.65],'LineWidth',1);
set(MUAPLabelAxes_Left,'XTick',[], ...
            'xTickLabel',[], ...
            'YLim',[0,(r+0.5)*(options.HeightSpacing+options.HeightPadding)+options.HeightPadding], ...
            'YTick',MUAP_Label_Tick_Values, ...
            'YTickLabel',MU_Labels, ...
            'FontSize',options.FontSize);
ylabel(MUAPLabelAxes_Left, options.YLabelLeft, 'FontSize', options.FontSize);

set(fig,'Color',options.BackColor, ...
      'WindowState','maximized', ...
      'PaperPositionMode','auto');
set(IDRAxes_Right,'Color','None', ...
    'XColor',options.FrontColor, ...
    'YColor',options.FrontColor, ...
    'ZColor',options.FrontColor, ...
    'FontName','Tahoma');

set(MUAPLabelAxes_Left,'Color','None', ...
    'XColor',options.FrontColor, ...
    'YColor',options.FrontColor, ...
    'ZColor',options.FrontColor, ...
    'FontName','Tahoma');

ckc.addTitles(MUAPLabelAxes_Left,options.Title,options.Subtitle);

end