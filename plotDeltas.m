function fig = plotDeltas(Delta)
%PLOTDELTAS  Plot arrays of nearest neighbor-intervals
%
% Syntax:
%   fig = ckc.plotDeltas(Delta);
%
% Inputs:
%   Delta - Returned by `ckc.compute_nearest_deltas`
%
% Output:
%   fig   - Figure handle with axes containing graphic of interval distribution.
%
% See also: Contents
fig = figure('Color', 'w', 'Name', 'Nearest-Neighbor Deltas');
ax = axes(fig,'NextPlot','add','FontName','Tahoma');

k = size(Delta,1);
title(ax,sprintf('\\Delta(%d) (N = %d)', k, size(Delta,2)), ...
    "FontName",'Tahoma','Color','k');
x = 1:size(Delta,2);
x_offset = linspace(0, 0.5, k);
cdata = colormap(jet(k));
for ii = 1:k
    line(ax, x + x_offset(ii), Delta(ii,:), 'LineStyle', 'none', 'Marker', 'o', ...
        'MarkerFaceColor', cdata(ii,:), 'MarkerEdgeColor', cdata(ii,:), ...
        'DisplayName', sprintf('k = %d', ii));

end
legend(ax,'FontName','Tahoma','TextColor','black');
xlabel(ax, 'MUPulse Number', 'FontName','Tahoma','Color','k');
ylabel(ax,'Samples to k-th Neighbor', 'FontName','Tahoma','Color','k');

ymed = median(Delta(:));
ymed_dev = median(abs(Delta(:) - ymed));
ylim(ax, [ymed-3*ymed_dev, ymed+3*ymed_dev]);

yline(ax, median(Delta(1,:)), 'Color', cdata(1,:), 'LineStyle', ':', 'Label', sprintf('N = %d', round(median(Delta(1,:)))), ...
    'DisplayName', 'Median \Delta');
ax.XLim(2) = ax.XLim(2)*1.25;

end