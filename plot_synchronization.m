function [fig, nUnits, goodCombinations, S, bins, N] = plot_synchronization(data, n, options)
%PLOT_SYNCHRONIZATION Plot the synchronization between MUAP pairs.
%
% Example:
%   [data,metadata,n] = io.load_cleaned_decomposition(14);
%   [fig, nUnits, goodCombinations, S, bins, N] = ckc.plot_synchronization(data,n);

arguments
    data
    n
    options.FigureParameters = {};
    options.MUAPOffset = 0;
    options.CMapData = [];
    options.MinNMuaps = 20;
    options.XScalar = 0.5; % For 2kHz to convert to milliseconds
    options.HorizontalUnits = 'Bin (ms)';
    options.VerticalUnits = 'Synchronization (a.u.)';
    options.Title = "";
    options.Subtitle = "";
    options.SynchronizationThreshold = 0.75;
end
[S, bins, N] = ckc.compute_synchronization(data);
goodCombinations = n.Pulses(:,[1,2]);
fig = figure('Color','w',...
    'Name','Synchronization MUAPs',...
    'WindowState','maximized', ...
    options.FigureParameters{:}); 
if isempty(options.CMapData)
    cmapData = jet(size(S,2));
else
    cmapData = options.CMapData;
end

bins = bins * options.XScalar;
dx = mean(diff(bins));
barCenters = bins(2:end)-0.5*dx;
val = min(abs(barCenters));
idx = abs(barCenters)==val; % In case it's non-unique


L = tiledlayout(fig, size(S,1)-1,size(S,2),'TileIndexing','columnmajor');
iLabel = (1:size(S,2)) + options.MUAPOffset;

xl = [bins(1), bins(end)];
ax = gobjects(size(S,1)-1,size(S,2));
nUnits = 0;
iRemove = [false(size(goodCombinations,1)-1,1);true];
for ii = 1:(size(S,1)-1)
    for ik = ii:size(S,2)
        iAx = sub2ind([size(S,1)-1,size(S,2)],ii,ik);
        ax(ii,ik) = nexttile(L,iAx,[1 1]);
        set(ax(ii,ik),'NextPlot','add','FontName','Tahoma','YLim',[0 1]);
        if ii == ik
            set(ax(ii,ik),'XLim',[0 1],'XColor','none','YColor','none');
            text(ax(ii,ik),0,0.5,sprintf('G-%d | U-%d',n.Pulses(ii,1),n.Pulses(ii,2)), ...
                'FontName','Tahoma','Color',cmapData(ii,:));
            text(ax(ii,ik),0,-1.5,sprintf('(N = %d)', N(ii)), ...
                'FontName','Tahoma','Color','k');
            isDuplicate = sum(iLabel==iLabel(ii))>1;
            tooFew = N(ii) < options.MinNMuaps;
            if isDuplicate || tooFew
                line(ax(ii,ik),[0 1],[0.5 0.5],'Color','r','LineWidth',1.5);
                set(ax(ii,ik),'Color','k');
                for ikk = 1:(ii-1)
                    line(ax(ikk,ii),xl,[0 1],'Color','r','LineWidth',1);
                    line(ax(ikk,ii),xl,[1 0],'Color','r','LineWidth',1);
                    set(ax(ikk,ii),'Color','k');
                end
                iRemove(ii) = true;
            else
                nUnits = nUnits + 1;
            end
        else
            set(ax(ii,ik),'XLim',xl);
            if any(S{ii,ik}(idx) > options.SynchronizationThreshold)
                cmapData(ik,:) = cmapData(ii,:);
                iLabel(ik) = iLabel(ii);
            end
            bar(ax(ii,ik),barCenters,S{ii,ik},1,...
                'EdgeColor','none','FaceColor',cmapData(ik,:)); 
            if isDuplicate || tooFew
                line(ax(ii,ik),xl,[0 1],'Color','r','LineWidth',1);
                line(ax(ii,ik),xl,[1 0],'Color','r','LineWidth',1);
                if tooFew
                    set(ax(ii,ik),'Color',[0.6 0.6 0.6]);
                else
                    set(ax(ii,ik),'Color','k');
                end
            end
            
        end
    end
end

xlabel(L, sprintf('Target | %s', options.HorizontalUnits), 'FontName','Tahoma');
ylabel(L, sprintf('Trigger | %s', options.VerticalUnits), 'FontName','Tahoma');
if strlength(options.Title) > 0
    title(L, options.Title, 'FontName','Tahoma','Color','k');
end
if strlength(options.Subtitle) > 0
    subtitle(L, options.Subtitle, 'FontName','Tahoma','Color',[0.65 0.65 0.65]);
else
    subtitle(L, sprintf('Sync Threshold: %4.2f | Unique MUAPs = %d', options.SynchronizationThreshold,nUnits), 'FontName','Tahoma','Color',[0.65 0.65 0.65]);
end
goodCombinations(iRemove,:) = [];
end