% Resolution Comparison — Max Effective Stress/Strain
% .mat files produced from postProcessing.mat script

clc
clear

d010 = load('F:\NCUR Simulations\4.2_0.1Res_NCPlaque\maxCalcResults_res010.mat');
d012 = load('F:\NCUR Simulations\8.4_0.12Res\maxCalcResults_res012.mat');

key010 = [d010.kappa(:), d010.rotation(:), d010.stenosis(:)];
key012 = [d012.kappa(:), d012.rotation(:), d012.stenosis(:)];

tissueFields = {'total', 'plaque', 'media', 'adventitia'};
tissueLabels = {'Total', 'Plaque', 'Media', 'Adventitia'};
quantFields  = {'effectiveStrain', 'effectiveStress'};
quantLabels  = {'EffectiveStrain', 'EffectiveStress'};

nMatched = 0;
pctDiff  = struct();
for qi = 1:length(quantFields)
    for ti = 1:length(tissueFields)
        pctDiff.(quantFields{qi}).(tissueFields{ti}) = [];
    end
end
matchedParams = [];

for i = 1:size(key010,1)
    match = find(all(abs(key012 - key010(i,:)) < 1e-6, 2), 1);
    if isempty(match)
        fprintf('No res0.12 match for kappa=%.0f rotation=%.0f stenosis=%.0f%% — skipping\n', ...
            key010(i,1), key010(i,2), key010(i,3)*100);
        continue
    end

    nMatched = nMatched + 1;
    matchedParams(nMatched,:) = key010(i,:);

    for qi = 1:length(quantFields)
        for ti = 1:length(tissueFields)
            qf = quantFields{qi};
            tf = tissueFields{ti};

            v1 = d010.sumData.(qf).(tf)(i);
            v2 = d012.sumData.(qf).(tf)(match);

            pct = abs((v2 - v1) / v1 * 100);
            pctDiff.(qf).(tf)(nMatched,1) = pct;
        end
    end
end

for qi = 1:length(quantFields)
    fprintf('-- %s --\n', quantLabels{qi});
    for ti = 1:length(tissueFields)
        vals = pctDiff.(quantFields{qi}).(tissueFields{ti});
        fprintf('  %s: Mean = %.2f%%, Max = %.2f%%\n', ...
            tissueLabels{ti}, mean(vals,'omitnan'), max(vals));
    end
    fprintf('\n');
end

save(fullfile(fileparts(mfilename('fullpath')), 'resolutionComparison.mat'), ...
     'matchedParams', 'pctDiff');

 
% Bar Charts (Max comparison)
resLabel010 = '0.10 mm';
resLabel012 = '0.12 mm';
color010 = [0.68 0.82 0.94];  % blue
color012 = [0.66 0.85 0.64];  % green
 
tissueOrder      = {'plaque','media','adventitia'};
tissueLabelsPlot = {'Plaque','Media','Adventitia'};
 
plotDefs = { ...
    struct('field','effectiveStress', 'title','Effective Stress', 'xlabel','Effective Stress (MPa)'), ...
    struct('field','effectiveStrain', 'title','Effective Strain', 'xlabel','Effective Strain') };
 
figure('Position',[100 100 1200 450], 'Color','white');
sgtitle('Mesh Resolution Comparison', ...
    'FontWeight','bold','FontSize',15);
 
for pd = 1:numel(plotDefs)
    def = plotDefs{pd};
    qf  = def.field;
 
    v010 = zeros(1,numel(tissueOrder));
    v012 = zeros(1,numel(tissueOrder));
    pctLabel = zeros(1,numel(tissueOrder));
 
    for ti = 1:numel(tissueOrder)
        tf = tissueOrder{ti};
        vals = pctDiff.(qf).(tf);
        [maxPct, idx] = max(vals);
 
        params = matchedParams(idx,:);
        i010 = find(all(abs(key010 - params) < 1e-6, 2), 1);
        i012 = find(all(abs(key012 - params) < 1e-6, 2), 1);
 
        v010(ti) = d010.sumData.(qf).(tf)(i010);
        v012(ti) = d012.sumData.(qf).(tf)(i012);
        pctLabel(ti) = maxPct;
    end
 
    ax = subplot(1, 2, pd);
    hold(ax, 'on');
 
    y = 1:numel(tissueOrder);
    barHeight = 0.32;
    b012 = barh(ax, y + barHeight/2, v012, barHeight, 'FaceColor', color012, 'EdgeColor','k','FaceAlpha',0.75);
    b010 = barh(ax, y - barHeight/2, v010, barHeight, 'FaceColor', color010, 'EdgeColor','k','FaceAlpha',0.75);

    xRange = max([v010 v012]);
    labelOffset = 0.03 * xRange;
 
    for ti = 1:numel(tissueOrder)
        xEnd = max(v010(ti), v012(ti));
        text(ax, xEnd + labelOffset, y(ti), sprintf('%.2f%%', pctLabel(ti)), ...
            'FontSize',11);
    end
 
    set(ax, 'YTick', y, 'YTickLabel', tissueLabelsPlot, 'FontSize',11, ...
        'TickDir','out','Box','off');
    xlabel(ax, def.xlabel, 'FontSize',11);
    title(ax, def.title, 'FontWeight','bold','FontSize',13);
    if pd == 2
        lgd = legend(ax, [b010 b012], {resLabel010, resLabel012}, ...
            'Orientation', 'vertical', 'Box', 'off', 'FontSize', 11);
        lgd.Units = 'normalized';
        lgd.Position(1) = ax.Position(1) + ax.Position(3) - lgd.Position(3) + 0.01;
        lgd.Position(2) = ax.Position(2) + ax.Position(4) - 0.035;
    end
    grid on;
end
