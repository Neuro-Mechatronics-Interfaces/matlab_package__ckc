function fig = plot_muap_raster(data, iGrid, iMUAP, options)
%PLOT_MUAP_RASTER Plots the raster for specified grid/muap combos using cleaned CKC data.

arguments
    data
    iGrid
    iMUAP
    options.SampleRate = 2000;
    options.GridID (1,:) string = ["DistExt"; "ProxExt"; "DistFlex"; "ProxFlex"];
    options.GridColor (:,3) = validatecolor(["#E9502C"; "#3A71E7"; "#F929A1"; "#A7F900"],"multiple");
    options.RefSignalIndex (1,2) {mustBePositive, mustBeInteger} = [1, 2]; % First is grid, second is channel.
    options.RefEMGIndex (1,2) {mustBePositive, mustBeInteger} = [1, 19]; 
    options.Title {mustBeTextScalar} = "";
    options.Subtitle {mustBeTextScalar} = "";
end

fig = figure('Color', 'w', 'Name', 'CKC Demuse Cleaned Raster');
switch getenv("COMPUTERNAME")
    case "MAX_LENOVO"
        fig.Position = [488   127   560   631];
end
ax = axes(fig,'NextPlot','add','FontName','Tahoma','YLim',[0,numel(iMUAP)+2]);
iFirst = nan(size(iMUAP));
pulses = cell(size(iMUAP));
for ii = 1:numel(iGrid)
    iFirst(ii) = min(data(iGrid(ii)).MUPulses{iMUAP(ii)});
    pulses{ii} = data(iGrid(ii)).MUPulses{iMUAP(ii)};
end
[~,idx] = sort(iFirst,'ascend');
pulses = pulses(idx);
iGrid = iGrid(idx);
iMUAP = iMUAP(idx);
lab = strings(numel(iGrid),1);
for ii = 1:numel(iGrid)
    lab(ii) = sprintf("%s-%d",options.GridID(iGrid(ii)),iMUAP(ii));
end
ref = data(options.RefSignalIndex(1)).ref_signal(options.RefSignalIndex(2),:);
t = 0:(1/options.SampleRate):((numel(ref)-1)/options.SampleRate);
set(ax,'XLim',[t(1),t(end)],'YTick',[],'YColor','none');
for ii = 1:numel(pulses)
    x = t(pulses{ii});
    xx = [x; x; nan(1,numel(x))];
    yy = [ones(2,numel(x)).*[ii-0.9; ii-0.1]; nan(1,numel(x))];
    line(ax,xx(:),yy(:),'LineWidth',1.25,'Color',options.GridColor(iGrid(ii),:));
    text(ax,0,ii-0.5,lab(ii),'FontName','Tahoma','Color',options.GridColor(iGrid(ii),:));
end
plot(ax,t,ref./max(ref)+numel(pulses),'Color',[0.65 0.65 0.65],'LineWidth',1.25);
data_ref = data(options.RefEMGIndex(1)).SIG{options.RefEMGIndex(2)};
[b,a] = butter(3,100./(options.SampleRate/2),'high');
data_ref = filtfilt(b,a,data_ref);
nMin = min(numel(t),numel(data_ref));
plot(ax,t(1:nMin),data_ref(1:nMin)./(2.*max(abs(data_ref(1:nMin)))) + numel(pulses) + 1.5,'Color',options.GridColor(options.RefEMGIndex(1),:),'LineWidth',1.0);
xlabel(ax,'Time (s)', 'FontName','Tahoma','Color','k');
text(ax,0,numel(pulses)+1.8,sprintf('%s: Channel-%d', options.GridID(options.RefEMGIndex(1)), options.RefEMGIndex(2)), ...
    'FontName','Tahoma','Color',options.GridColor(options.RefEMGIndex(1),:));
if strlength(options.Title) > 0
    title(ax,options.Title,'FontName','Tahoma','Color','k');
end

if strlength(options.Subtitle) > 0
    subtitle(ax,options.Subtitle,'FontName','Tahoma','Color',[0.65 0.65 0.65]);
end
end