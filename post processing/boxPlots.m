%% Load data
close all;
clear;
clc;

% Paths and Data
path = "F:\NCUR Simulations\8.4_0.12Res";
plottype = 'Max';
dataDirection = 'Axial';
summaryMat = fullfile(path, [lower(plottype) 'CalcResults_res012.mat']);

R = load(summaryMat, ...
    'kappa', 'rotation', 'stenosis', 'resolution', 'edge_length', ...
    'sumData');

% Stress & Strain Calculations
R.rotation = wrapTo360(R.rotation);
for i=1:length(R.kappa)
    if R.kappa(i) < 0
        switch R.rotation(i)
            case 0
                R.rotation(i) = 180;
        end
    end
end

kappa    = unique(abs(R.kappa));
rot      = unique(wrapTo360(R.rotation));
stenosis = unique(R.stenosis)*100;

nK = length(kappa);
nR = length(rot);
nS = length(stenosis);

% Preallocate matrices
totalTensileStress      = zeros(nK, nR, nS);
plaqueTensileStress     = zeros(nK, nR, nS);
mediaTensileStress      = zeros(nK, nR, nS);
advenTensileStress      = zeros(nK, nR, nS);

totalCompressiveStress  = zeros(nK, nR, nS);
plaqueCompressiveStress = zeros(nK, nR, nS);
mediaCompressiveStress  = zeros(nK, nR, nS);
advenCompressiveStress  = zeros(nK, nR, nS);

totalTensileStrain      = zeros(nK, nR, nS);
plaqueTensileStrain     = zeros(nK, nR, nS);
mediaTensileStrain      = zeros(nK, nR, nS);
advenTensileStrain      = zeros(nK, nR, nS);

totalCompressiveStrain  = zeros(nK, nR, nS);
plaqueCompressiveStrain = zeros(nK, nR, nS);
mediaCompressiveStrain  = zeros(nK, nR, nS);
advenCompressiveStrain  = zeros(nK, nR, nS);

totalEffectiveStrain    = zeros(nK, nR, nS);
plaqueEffectiveStrain   = zeros(nK, nR, nS);
mediaEffectiveStrain    = zeros(nK, nR, nS);
advenEffectiveStrain    = zeros(nK, nR, nS);
 
totalEffectiveStress    = zeros(nK, nR, nS);
plaqueEffectiveStress   = zeros(nK, nR, nS);
mediaEffectiveStress    = zeros(nK, nR, nS);
advenEffectiveStress    = zeros(nK, nR, nS);

% Fill matrices
for i = 1:nK
    for j = 1:nR
        for k = 1:nS
            idx = (abs(R.kappa) == kappa(i)) & ...
                  (wrapTo360(R.rotation) == rot(j)) & ...
                  (R.stenosis*100 == stenosis(k));
            if any(idx)
                f = find(idx,1,'first');
                plaqueEffectiveStress(i,j,k)   = R.sumData.effectiveStress.plaque(f);
                mediaEffectiveStress(i,j,k)   = R.sumData.effectiveStress.media(f);
                advenEffectiveStress(i,j,k)   = R.sumData.effectiveStress.adventitia(f);
                totalEffectiveStress(i,j,k)   = R.sumData.effectiveStress.total(f);

                plaqueEffectiveStrain(i,j,k)   = R.sumData.effectiveStrain.plaque(f);
                mediaEffectiveStrain(i,j,k)   = R.sumData.effectiveStrain.media(f);
                advenEffectiveStrain(i,j,k)   = R.sumData.effectiveStrain.adventitia(f);
                totalEffectiveStrain(i,j,k)   = R.sumData.effectiveStrain.total(f);

                plaqueTensileStress(i,j,k)     = R.sumData.tensileStress.plaque(f);
                mediaTensileStress(i,j,k)      = R.sumData.tensileStress.media(f);
                advenTensileStress(i,j,k)      = R.sumData.tensileStress.adventitia(f);
                totalTensileStress(i,j,k)      = R.sumData.tensileStress.total(f);

                plaqueCompressiveStress(i,j,k) = R.sumData.compressiveStress.plaque(f);
                mediaCompressiveStress(i,j,k)  = R.sumData.compressiveStress.media(f);
                advenCompressiveStress(i,j,k)  = R.sumData.compressiveStress.adventitia(f);
                totalCompressiveStress(i,j,k)  = R.sumData.compressiveStress.total(f);

                plaqueTensileStrain(i,j,k)     = R.sumData.tensileStrain.plaque(f);
                mediaTensileStrain(i,j,k)      = R.sumData.tensileStrain.media(f);
                advenTensileStrain(i,j,k)      = R.sumData.tensileStrain.adventitia(f);
                totalTensileStrain(i,j,k)      = R.sumData.tensileStrain.total(f);

                plaqueCompressiveStrain(i,j,k) = R.sumData.compressiveStrain.plaque(f);
                mediaCompressiveStrain(i,j,k)  = R.sumData.compressiveStrain.media(f);
                advenCompressiveStrain(i,j,k)  = R.sumData.compressiveStrain.adventitia(f);
                totalCompressiveStrain(i,j,k)  = R.sumData.compressiveStrain.total(f);
            end
        end
    end
end

% Colors
cPlaque = [0.27 0.52 0.90];   % blue
cMedia  = [0.18 0.65 0.18];   % green
cAdven  = [0.90 0.45 0.18];   % orange

%% Generate 6 separate figures: 3 parameters x 2 metrics

paramNames  = {'Stenosis Severity', 'Curvature Direction', 'Curvature Magnitude'};
paramYLabel = {'Stenosis Severity (%)', 'Curvature Direction (°)', 'Curvature Magnitude (\kappa)'};
 
% Stress Box Plots
for p = 1:3
    makeSingleFigure(plaqueTensileStress, mediaTensileStress, advenTensileStress, ...
                      plaqueCompressiveStress, mediaCompressiveStress, advenCompressiveStress, ...
                      kappa, rot, stenosis, nK, nR, nS, cPlaque, cMedia, cAdven, ...
                      'stress', dataDirection, p, paramNames{p}, paramYLabel{p}, plottype);
end

% Strain Box Plots
for p = 1:3
    makeSingleFigure(plaqueTensileStrain, mediaTensileStrain, advenTensileStrain, ...
                      plaqueCompressiveStrain, mediaCompressiveStrain, advenCompressiveStrain, ...
                      kappa, rot, stenosis, nK, nR, nS, cPlaque, cMedia, cAdven, ...
                      'strain', dataDirection, p, paramNames{p}, paramYLabel{p}, plottype);
end

% Effective Stress Box Plots
for p = 1:3
    makeSingleEffectiveFigure(plaqueEffectiveStress, mediaEffectiveStress, advenEffectiveStress, ...
                      kappa, rot, stenosis, nK, nR, nS, cPlaque, cMedia, cAdven, ...
                      'Effective Stress', dataDirection, p, paramNames{p}, paramYLabel{p}, plottype);
end

% Effective Strain Box Plots
for p = 1:3
    makeSingleEffectiveFigure(plaqueEffectiveStrain, mediaEffectiveStrain, advenEffectiveStrain, ...
                      kappa, rot, stenosis, nK, nR, nS, cPlaque, cMedia, cAdven, ...
                      'Effective Strain', dataDirection, p, paramNames{p}, paramYLabel{p}, plottype);
end

%% Helper functions

function makeSingleFigure(pTen, mTen, aTen, pCom, mCom, aCom, ...
                    kappa, rot, stenosis, nK, nR, nS, cP, cM, cA, ...
                    metricType, dataDirection, paramIdx, paramName, paramYLabel, plottype)
 
    figure('Position', [100 100 1100 500], 'Color', 'white');
 
    allVals = [pTen(:); mTen(:); aTen(:); pCom(:); mCom(:); aCom(:)];
    allVals = allVals(allVals ~= 0 & ~isnan(allVals));
    maxVal  = max(abs(allVals));
    if maxVal == 0, maxVal = 1; end
    exponent = floor(log10(maxVal));
    scale    = 10^(-exponent);
    expLabel = num2str(-exponent);
 
    if strcmp(metricType, 'stress')
        xLabel = sprintf('Stress (\\times10^{-%s} MPa)', expLabel);
    else
        xLabel = sprintf('Strain (\\times10^{-%s})', expLabel);
    end
 
    allTen = [pTen(:); mTen(:); aTen(:)];
    allCom = [pCom(:); mCom(:); aCom(:)];
    allTen = allTen(allTen > 0 & ~isnan(allTen)) * scale;
    allCom = allCom(allCom < 0 & ~isnan(allCom)) * scale;
 
    pad   = 0.08;
    xMinC = min(allCom) * (1 + pad);
    xMaxC = max(allCom) * (1 - pad);
    xMinT = min(allTen) * (1 - pad);
    xMaxT = max(allTen) * (1 + pad);
 
    % Layout: single row, two panels
    leftMar  = 0.14;
    rightMar = 0.03;
    midGap   = 0.04;
    botMar   = 0.14;
    topMar   = 0.12;
    rowH     = 1 - botMar - topMar;
    yb       = botMar;
 
    totalFrac = 1 - leftMar - rightMar - midGap;
    compFrac  = totalFrac / 2;
    tenFrac   = totalFrac / 2;
 
    % Reshape data depending on which parameter is being plotted
    nObs  = nK * nR;
    nObs2 = nK * nS;
    nObs3 = nR * nS;
 
    switch paramIdx
        case 1  % vs stenosis severity
            dp_T = reshape(pTen, nObs, nS) * scale;
            dm_T = reshape(mTen, nObs, nS) * scale;
            da_T = reshape(aTen, nObs, nS) * scale;
            dp_C = reshape(pCom, nObs, nS) * scale;
            dm_C = reshape(mCom, nObs, nS) * scale;
            da_C = reshape(aCom, nObs, nS) * scale;
            groupLabels = arrayfun(@(s) sprintf('%d%%', s), stenosis, 'UniformOutput', false);
        case 2  % vs curvature direction (rotation)
            dp_T = reshape(permute(pTen,[1 3 2]), nObs2, nR) * scale;
            dm_T = reshape(permute(mTen,[1 3 2]), nObs2, nR) * scale;
            da_T = reshape(permute(aTen,[1 3 2]), nObs2, nR) * scale;
            dp_C = reshape(permute(pCom,[1 3 2]), nObs2, nR) * scale;
            dm_C = reshape(permute(mCom,[1 3 2]), nObs2, nR) * scale;
            da_C = reshape(permute(aCom,[1 3 2]), nObs2, nR) * scale;
            groupLabels = arrayfun(@(r) sprintf('%d°', r), rot, 'UniformOutput', false);
        case 3  % vs curvature magnitude (kappa)
            dp_T = reshape(permute(pTen,[2 3 1]), nObs3, nK) * scale;
            dm_T = reshape(permute(mTen,[2 3 1]), nObs3, nK) * scale;
            da_T = reshape(permute(aTen,[2 3 1]), nObs3, nK) * scale;
            dp_C = reshape(permute(pCom,[2 3 1]), nObs3, nK) * scale;
            dm_C = reshape(permute(mCom,[2 3 1]), nObs3, nK) * scale;
            da_C = reshape(permute(aCom,[2 3 1]), nObs3, nK) * scale;
            groupLabels = arrayfun(@(k) sprintf('%d', k), kappa, 'UniformOutput', false);
    end
 
    % Compressive panel
    axC = axes('Position', [leftMar, yb, compFrac, rowH]);
    doGroupedBoxplot(axC, dp_C, dm_C, da_C, groupLabels, cP, cM, cA);
    delete(findobj(axC, 'Type', 'ConstantLine'));
    xlim(axC, [xMinC, xMaxC]);
    set(axC, 'XDir', 'reverse', 'Box', 'off');
    xlabel(axC, xLabel, 'FontSize', 11);
    ylabel(axC, paramYLabel, 'FontSize', 11);
    title(axC, 'Compressive', 'FontSize', 11, 'FontWeight', 'bold')
 
    % Tensile panel
    axT = axes('Position', [leftMar + compFrac + midGap, yb, tenFrac, rowH]);
    h = doGroupedBoxplot(axT, dp_T, dm_T, da_T, groupLabels, cP, cM, cA);
    delete(findobj(axT, 'Type', 'ConstantLine'));
    xlim(axT, [xMinT, xMaxT]);
    set(axT, 'YTickLabel', [], 'Box', 'off');
    xlabel(axT, xLabel, 'FontSize', 11);
    title(axT, 'Tensile', 'FontSize', 11, 'FontWeight', 'bold')
 
    % Divider line between panels
    xDiv = leftMar + compFrac + midGap/2;
    annotation('line', [xDiv xDiv], [yb, yb + rowH], ...
        'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);
 
    % Title and legend
    metricLabel = upper(metricType(1)) + string(metricType(2:end)); % 'Stress' or 'Strain'
    fullLabel   = sprintf('%s %s %s', plottype, dataDirection, metricLabel);  % e.g. 'Max Axial Stress'
    sgtitle(sprintf('%s — %s', fullLabel, paramName), 'FontSize', 13, 'FontWeight', 'bold');
    
    if paramIdx == 1
        lgd = legend(axT, h, {'Plaque', 'Media', 'Adventitia'}, ...
            'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9);
        lgd.Units = 'normalized';
        lgd.Position(1) = 1 - rightMar - lgd.Position(3);
        lgd.Position(2) = yb + rowH + 0.04;
    end
end

%% makeSingleEffectiveFigure
function makeSingleEffectiveFigure(pEff, mEff, aEff, ...
                    kappa, rot, stenosis, nK, nR, nS, cP, cM, cA, ...
                    metricLabel, dataDirection, paramIdx, paramName, paramYLabel, plottype)
 
    figure('Position', [100 100 1100 500], 'Color', 'white');
 
    allVals = [pEff(:); mEff(:); aEff(:)];
    allVals = allVals(allVals ~= 0 & ~isnan(allVals));
    maxVal  = max(abs(allVals));
    if maxVal == 0, maxVal = 1; end
    exponent = floor(log10(maxVal));
    scale    = 10^(-exponent);
    expLabel = num2str(-exponent);
 
    if contains(lower(metricLabel), 'stress')
        xLabel = sprintf('%s (\\times10^{-%s} MPa)', metricLabel, expLabel);
    else
        xLabel = sprintf('%s (\\times10^{-%s})', metricLabel, expLabel);
    end
 
    allData = allVals * scale;
    pad   = 0.08;
    xMinV = 0;
    xMaxV = max(allData) * (1 + pad);
 
    % Layout: single panel
    leftMar  = 0.20;
    rightMar = 0.06;
    botMar   = 0.14;
    topMar   = 0.12;
    rowH     = 1 - botMar - topMar;
    colW     = 1 - leftMar - rightMar;
    yb       = botMar;
 
    % Reshape data depending on which parameter is being plotted
    nObs  = nK * nR;
    nObs2 = nK * nS;
    nObs3 = nR * nS;
 
    switch paramIdx
        case 1  % vs stenosis severity
            dp = reshape(pEff, nObs, nS) * scale;
            dm = reshape(mEff, nObs, nS) * scale;
            da = reshape(aEff, nObs, nS) * scale;
            groupLabels = arrayfun(@(s) sprintf('%d%%', s), stenosis, 'UniformOutput', false);
        case 2  % vs curvature direction (rotation)
            dp = reshape(permute(pEff,[1 3 2]), nObs2, nR) * scale;
            dm = reshape(permute(mEff,[1 3 2]), nObs2, nR) * scale;
            da = reshape(permute(aEff,[1 3 2]), nObs2, nR) * scale;
            groupLabels = arrayfun(@(r) sprintf('%d°', r), rot, 'UniformOutput', false);
        case 3  % vs curvature magnitude (kappa)
            dp = reshape(permute(pEff,[2 3 1]), nObs3, nK) * scale;
            dm = reshape(permute(mEff,[2 3 1]), nObs3, nK) * scale;
            da = reshape(permute(aEff,[2 3 1]), nObs3, nK) * scale;
            groupLabels = arrayfun(@(k) sprintf('%d', k), kappa, 'UniformOutput', false);
    end
 
    ax = axes('Position', [leftMar, yb, colW, rowH]);
    h = doGroupedBoxplot(ax, dp, dm, da, groupLabels, cP, cM, cA);
    delete(findobj(ax, 'Type', 'ConstantLine'));
    xlim(ax, [xMinV, xMaxV]);
    set(ax, 'Box', 'off');
    xlabel(ax, xLabel, 'FontSize', 11);
    ylabel(ax, paramYLabel, 'FontSize', 11);
 
    fullLabel = sprintf('%s %s %s', plottype, dataDirection, metricLabel);  % e.g. 'Max Axial Effective Stress'
    sgtitle(sprintf('%s — %s', fullLabel, paramName), 'FontSize', 13, 'FontWeight', 'bold');
 
    if paramIdx == 1
        lgd = legend(ax, h, {'Plaque', 'Media', 'Adventitia'}, ...
            'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9);
        lgd.Units = 'normalized';
        lgd.Position(1) = 1 - rightMar - lgd.Position(3);
        lgd.Position(2) = yb + rowH + 0.04;
    end
end
 
%% doGroupedBoxplot
function h = doGroupedBoxplot(ax, dp, dm, da, groupLabels, cP, cM, cA)
    axes(ax); hold on;
    nGroups  = size(dp, 2);
    boxWidth = 0.18;
    offsets  = [-1, 0, 1] * boxWidth;
    colors   = {cP, cM, cA};
    data     = {dp, dm, da};
 
    for s = 1:3
        for g = 1:nGroups
            yc  = g + offsets(s);
            col = data{s}(:,g);
            col = col(col ~= 0 & ~isnan(col));
 
            if isempty(col), continue; end
 
            q = quantile(col, [0.0 0.25 0.5 0.75 1.0]);
 
            % Horizontal box
            yBox = yc + boxWidth/2 * [-1 -1  1  1 -1];
            xBox = [q(2) q(4) q(4) q(2) q(2)];
            patch(xBox, yBox, colors{s}, ...
                'FaceAlpha', 0.4, 'EdgeColor', colors{s}, 'LineWidth', 1.2);
 
            % Median line
            line([q(3) q(3)], [yc - boxWidth/2, yc + boxWidth/2], ...
                'Color', 'k', 'LineWidth', 1.8);
 
            % Whiskers
            line([q(1) q(2)], [yc yc], 'Color', colors{s}, 'LineWidth', 1);
            line([q(4) q(5)], [yc yc], 'Color', colors{s}, 'LineWidth', 1);
 
            % Whisker caps
            capW = boxWidth * 0.3;
            line([q(1) q(1)], [yc-capW yc+capW], 'Color', colors{s}, 'LineWidth', 1);
            line([q(5) q(5)], [yc-capW yc+capW], 'Color', colors{s}, 'LineWidth', 1);
        end
    end
 
    % Legend handles
    h(1) = bar(nan, nan, 'FaceColor', cP, 'FaceAlpha', 0.3, 'EdgeColor', cP);
    h(2) = bar(nan, nan, 'FaceColor', cM, 'FaceAlpha', 0.3, 'EdgeColor', cM);
    h(3) = bar(nan, nan, 'FaceColor', cA, 'FaceAlpha', 0.3, 'EdgeColor', cA);
 
    set(ax, 'YTick', 1:nGroups, 'YTickLabel', groupLabels, ...
        'FontSize', 9, 'Box', 'off', 'TickDir', 'out');
    ylim([0.5, nGroups + 0.5]);
    xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    grid on;
end