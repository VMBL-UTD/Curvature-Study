%% Load data
close all;
clear;
clc;

% Paths and Data
path = "F:\NCUR Simulations\8.4_0.12Res";
plottype = 'Max';
datatype = 'Axial';
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

%% Axial Stress (Tensile + Compressive)
makeFigure(plaqueTensileStress, mediaTensileStress, advenTensileStress, ...
           plaqueCompressiveStress, mediaCompressiveStress, advenCompressiveStress, ...
           kappa, rot, stenosis, nK, nR, nS, cPlaque, cMedia, cAdven, 'stress');

%% Axial Strain (Tensile + Compressive)
makeFigure(plaqueTensileStrain, mediaTensileStrain, advenTensileStrain, ...
           plaqueCompressiveStrain, mediaCompressiveStrain, advenCompressiveStrain, ...
           kappa, rot, stenosis, nK, nR, nS, cPlaque, cMedia, cAdven, 'strain');


%% makeFigure
function makeFigure(pTen, mTen, aTen, pCom, mCom, aCom, ...
                    kappa, rot, stenosis, nK, nR, nS, cP, cM, cA, datatype)

    figure('Position', [100 100 1000 1200], 'Color', 'white');

    allVals = [pTen(:); mTen(:); aTen(:); pCom(:); mCom(:); aCom(:)];
    allVals = allVals(allVals ~= 0 & ~isnan(allVals));
    maxVal  = max(abs(allVals));
    if maxVal == 0, maxVal = 1; end
    exponent = floor(log10(maxVal));
    scale    = 10^(-exponent);
    expLabel = num2str(-exponent);

    if strcmp(datatype, 'stress')
        xLabel = sprintf('Stress (\\times10^{-%s} MPa)', expLabel);
    else
        xLabel = sprintf('Strain (\\times10^{-%s})', expLabel);
    end

    keepFL = @(x) x;

    % Data ranges
    allTen = [pTen(:); mTen(:); aTen(:)];
    allCom = [pCom(:); mCom(:); aCom(:)];
    allTen = allTen(allTen > 0 & ~isnan(allTen)) * scale;
    allCom = allCom(allCom < 0 & ~isnan(allCom)) * scale;

    pad   = 0.08;
    xMinC = min(allCom) * (1 + pad);
    xMaxC = max(allCom) * (1 - pad);
    xMinT = min(allTen) * (1 - pad);
    xMaxT = max(allTen) * (1 + pad);

    % Layout
    leftMar  = 0.12;
    rightMar = 0.03;
    midGap   = 0.04;
    rowH     = 0.26;
    rowGap   = 0.06;
    botMar   = 0.08;

    totalFrac = 1 - leftMar - rightMar - midGap;
    compFrac  = totalFrac / 2;
    tenFrac   = totalFrac / 2;

    rowBottoms = [botMar + 2*(rowH+rowGap), botMar + (rowH+rowGap), botMar];

    % Prepare reshaped data
    nObs  = nK * nR;
    nObs2 = nK * nS;
    nObs3 = nR * nS;

    allDP_T = { keepFL(reshape(pTen,                 nObs,  nS) * scale), ...
                keepFL(reshape(permute(pTen,[1 3 2]), nObs2, nR) * scale), ...
                keepFL(reshape(permute(pTen,[2 3 1]), nObs3, nK) * scale) };
    allDM_T = { keepFL(reshape(mTen,                 nObs,  nS) * scale), ...
                keepFL(reshape(permute(mTen,[1 3 2]), nObs2, nR) * scale), ...
                keepFL(reshape(permute(mTen,[2 3 1]), nObs3, nK) * scale) };
    allDA_T = { keepFL(reshape(aTen,                 nObs,  nS) * scale), ...
                keepFL(reshape(permute(aTen,[1 3 2]), nObs2, nR) * scale), ...
                keepFL(reshape(permute(aTen,[2 3 1]), nObs3, nK) * scale) };

    allDP_C = { keepFL(reshape(pCom,                 nObs,  nS) * scale), ...
                keepFL(reshape(permute(pCom,[1 3 2]), nObs2, nR) * scale), ...
                keepFL(reshape(permute(pCom,[2 3 1]), nObs3, nK) * scale) };
    allDM_C = { keepFL(reshape(mCom,                 nObs,  nS) * scale), ...
                keepFL(reshape(permute(mCom,[1 3 2]), nObs2, nR) * scale), ...
                keepFL(reshape(permute(mCom,[2 3 1]), nObs3, nK) * scale) };
    allDA_C = { keepFL(reshape(aCom,                 nObs,  nS) * scale), ...
                keepFL(reshape(permute(aCom,[1 3 2]), nObs2, nR) * scale), ...
                keepFL(reshape(permute(aCom,[2 3 1]), nObs3, nK) * scale) };

    allLabels = { arrayfun(@(s) sprintf('%d%%', s), stenosis, 'UniformOutput', false), ...
                  arrayfun(@(r) sprintf('%d°',  r), rot,      'UniformOutput', false), ...
                  arrayfun(@(k) sprintf('%d',   k), kappa,    'UniformOutput', false) };

    axC = gobjects(3,1);
    axT = gobjects(3,1);

    for row = 1:3
        yb = rowBottoms(row);

        % Compressive panel
        axC(row) = axes('Position', [leftMar, yb, compFrac, rowH]);
        doGroupedBoxplot(axC(row), allDP_C{row}, allDM_C{row}, allDA_C{row}, allLabels{row}, cP, cM, cA);
        delete(findobj(axC(row), 'Type', 'ConstantLine'));
        xlim(axC(row), [xMinC, xMaxC]);
        set(axC(row), 'XDir', 'reverse', 'Box', 'off');
        if row ~= 3
            set(axC(row), 'XTickLabel', []);
        end
        hold(axC(row), 'off');

        % Tensile panel
        axT(row) = axes('Position', [leftMar + compFrac + midGap, yb, tenFrac, rowH]);
        h = doGroupedBoxplot(axT(row), allDP_T{row}, allDM_T{row}, allDA_T{row}, allLabels{row}, cP, cM, cA);
        delete(findobj(axT(row), 'Type', 'ConstantLine'));
        xlim(axT(row), [xMinT, xMaxT]);
        set(axT(row), 'YTickLabel', [], 'Box', 'off');
        if row ~= 3
            set(axT(row), 'XTickLabel', []);
        end
        hold(axT(row), 'off');

        if row == 3
            xlabel(axC(row), xLabel, 'FontSize', 11);
            xlabel(axT(row), xLabel, 'FontSize', 11);
        end

        % Divider line between panels
        xDiv = leftMar + compFrac + midGap/2;
        annotation('line', [xDiv xDiv], [yb, yb + rowH], ...
            'Color', [0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth', 1);
    end

    linkaxes(axC, 'x');
    linkaxes(axT, 'x');

    legend(axT(1), h, {'Plaque', 'Media', 'Adventitia'}, 'Location', 'northeast');
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