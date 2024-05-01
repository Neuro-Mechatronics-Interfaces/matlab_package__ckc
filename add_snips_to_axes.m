function add_snips_to_axes(ax, snips, tsnip, options)
%ADD_SNIPS_TO_AXES  Add snips from cleaned CKC signals to axes to plot templates/examples
%
% Syntax:
%   ckc.add_snips_to_axes(ax,snips,tsnip,'Name',value,...);
%
% Inputs:
%   ax - The axes handle to add the plot onto
%   snips - The cell array of waveform snippets that are "extended samples"
%           around each MUAP pulse detected in CKC/DEMUSE.
%   tsnip - The time vector corresponding to `snips`
arguments
    ax
    snips
    tsnip
    options.AddScaleBar    (1,1) logical = true;
    options.FormatAxes     (1,1) logical = true;
    options.PlotMean       (1,1) logical = true;
    options.PlotIndividual (1,1) logical = true;
    options.NIndividualMax (1,1) {mustBePositive,mustBeInteger} = 25;
    options.NGrid (1,1) {mustBePositive, mustBeInteger} = 4;
    options.NumRows (1,1) {mustBePositive, mustBeInteger} = 16;
    options.NumColumns (1,1) {mustBePositive, mustBeInteger} = 16;
    options.NPerGrid (1,:) {mustBePositive, mustBeInteger} = [64, 64, 64, 64];
    options.GridColor (:,3) = validatecolor(["#E9502C"; "#3A71E7";  "#A7F900"; "#F929A1" ],"multiple");
    options.GridOffset (1,1) double = 1.5; % Times usual offset
    options.GridQuadrantMap (:,2) {mustBeInteger} = [0 0; 1 0; 0 1; 1 1];
    options.YScale (1,1) double = 50; % microvolts
    options.YUnits {mustBeTextScalar} = '\muV';
    options.XScalar (1,1) double = 1.1; % Times the total time (width) of each snippet
    options.XUnits {mustBeTextScalar} = 'ms';
    options.DebugMode (1,1) logical = false;
end

offsetX = options.XScalar * (tsnip(end)-tsnip(1));
offsetY = options.YScale;
gridSpacingX = options.GridOffset*offsetX + 8*offsetX;
gridSpacingY = options.GridOffset*offsetY + 8*offsetY;
nGrid = options.NGrid;
if nGrid > size(options.GridColor,1)
    error("Must have at least as many rows (%d) in GridColor (currently: %d) option as there are grids (cell elements in snips)!", nGrid, size(options.GridColor,1));
end
nRows = options.NumRows;
nCols = options.NumColumns;
xl = [tsnip(1)-offsetX, 2*gridSpacingX-offsetX];
yl = [-offsetY, 2*gridSpacingY-offsetY];
if options.FormatAxes
    set(ax,'NextPlot','add','XColor','none','YColor','none','XLim',xl,'YLim',yl);
    switch getenv("COMPUTERNAME")
        case 'MAX_LENOVO'
            set(get(ax,'Parent'),'Position',[100   164   789   643]);
    end
end

nChTotal = nRows * nCols;
iSample = randsample(size(snips,1),min(size(snips,1),options.NIndividualMax),false);
for iCh = 1:nChTotal
    % [iRow, iCol] = ind2sub([nRows, nCols], iCh);
    % iGrid = 2*floor((iRow-1)/8) + floor((iCol-1)/8) + 1;
    iGrid = floor((iCh-1)/64)+1;
    xGridOffset = gridSpacingX*options.GridQuadrantMap(iGrid,1);
    yGridOffset = gridSpacingY*options.GridQuadrantMap(iGrid,2);
    gridCh = rem(iCh-1,64)+1;
    xChannelOffset = offsetX*floor((gridCh-1)/8);
    yChannelOffset = offsetY*rem((gridCh-1),8);
    if options.DebugMode
        fprintf(1,'iCh = %d | xyGridOffset = <%f, %f> | xyChannelOffset = <%f, %f>\n', iCh, xGridOffset, yGridOffset, xChannelOffset, yChannelOffset);
    end
    if options.PlotIndividual
        plot(ax, ...
            tsnip + xChannelOffset + xGridOffset, ...
            squeeze(snips(iSample,:,iCh)) + yChannelOffset + yGridOffset,...
            'LineWidth',0.5,'Color',[0.65 0.65 0.65]);
    end
    if options.PlotMean
        plot(ax, ...
            tsnip + xChannelOffset + xGridOffset, ...
            squeeze(mean(snips(:,:,iCh),1)) + yChannelOffset + yGridOffset,...
            'LineWidth',2.5,'Color',options.GridColor(iGrid,:));
    end
end
if options.AddScaleBar
    line(ax, [-offsetX, 0], [-offsetY, -offsetY], 'Color', 'k', 'LineWidth', 1);
    line(ax, [-offsetX, -offsetX], [-offsetY, 0], 'Color', 'k', 'LineWidth', 1);
    text(ax, -1.25*offsetX, -0.5*offsetY, sprintf('%d %s', round(options.YScale), options.YUnits), 'FontName','Tahoma', 'Color', 'k', 'HorizontalAlignment', 'right');
    text(ax, -0.5*offsetX, -1.25*offsetY, sprintf('%d %s', round(offsetX), options.XUnits), 'FontName', 'Tahoma', 'Color', 'k');
end

end