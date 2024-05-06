function [fig, MU_ID, recruitmentOrder] = plotIPTsFast(t, MUPulses, IPTs, MU_ID, recruitmentOrder, options)
%PLOTIPTSFAST  Create fast tiled plot of IPTs
%
% Syntax:
%   [fig, MU_ID, recruitmentOrder] = ckc.plotIPTsFast(t, MUPulses, IPTs, MU_ID, recruitmentOrder);
%
% Inputs:
%   t - (1 x nSamples) Time vector
%   MUPulses - Pulse instants (1 x nMUAPS cell array)
%   IPTs - IPTs from DEMUSE (pulse trains)
%   MU_ID - Motor Unit IDs for each IPT (optional; generates if left as [])
%   recruitmentOrder - Order of recruitment so that color scheme matches
%                       other plots (like `ckc.plotIDR` with 'Recruitment' SortMethod option)
%
% Output:
%   fig   - Generated figure handle
%   MU_ID - (Useful if you left input empty)
%   recruitmentOrder - (Useful if you left input empty)
%
% See also: Contents

arguments
    t (1,:) {mustBeNumeric}
    MUPulses (1,:) cell
    IPTs (:,:) {mustBeNumeric}
    MU_ID = [];
    recruitmentOrder = [];
    options.LineWidth (1,1) double = 0.5;
    options.XLabelString {mustBeTextScalar} = "Time (s)";
    options.YLabelString {mustBeTextScalar} = "IPT (a.u.)";
    options.Title {mustBeTextScalar} = "";
    options.Subtitle {mustBeTextScalar} = "";
    options.TitleColor (1,3) {mustBeInRange(options.TitleColor,0,1)} = [0 0 0];
    options.SubtitleColor (1,3) {mustBeInRange(options.SubtitleColor,0,1)} = [0.65 0.65 0.65];
    options.FontName {mustBeTextScalar} = 'Tahoma';
end

if isempty(MU_ID)
    MU_ID = ckc.generate_muap_id(MUPulses);
end

if isempty(recruitmentOrder)
    recruitmentOrder = ckc.sortByRecruitmentOrder(MUPulses);
end

fig = figure('Color','w','Position', [150, 100, 600, 500]); 
L = tiledlayout(fig,'flow'); 
cmapdata = jet(size(IPTs,1)); 
for ii = 1:size(IPTs,1) 
    ax = nexttile(L); 
    set(ax,'NextPlot','add', ...
        'FontName', options.FontName, 'FontSize', 12); 
    plot(ax, t, IPTs(recruitmentOrder(ii),:), ...
        'LineWidth', options.LineWidth, ...
        'Color', cmapdata(ii,:)); 
    title(ax, MU_ID{recruitmentOrder(ii)}, ...
        'FontName', options.FontName, ...
        'Color', cmapdata(ii,:)); 
end
xlabel(L, options.XLabelString, 'FontName', options.FontName, 'FontSize', 14);
ylabel(L, options.YLabelString, 'FontName', options.FontName, 'FontSize', 14);
ckc.addTitles(L,options.Title,options.Subtitle,options.TitleColor,options.SubtitleColor,'FontName',options.FontName);
end