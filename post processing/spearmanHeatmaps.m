% Spearman Rho Correlation Heatmaps
% Tensile and compressive stress correlation
 
% Paths and Loading Data
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
 
% Preallocate
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
                mediaEffectiveStress(i,j,k)    = R.sumData.effectiveStress.media(f);
                advenEffectiveStress(i,j,k)    = R.sumData.effectiveStress.adventitia(f);
                totalEffectiveStress(i,j,k)    = R.sumData.effectiveStress.total(f);
 
                plaqueTensileStress(i,j,k)     = R.sumData.tensileStress.plaque(f);
                mediaTensileStress(i,j,k)      = R.sumData.tensileStress.media(f);
                advenTensileStress(i,j,k)      = R.sumData.tensileStress.adventitia(f);
                totalTensileStress(i,j,k)      = R.sumData.tensileStress.total(f);
 
                plaqueCompressiveStress(i,j,k) = R.sumData.compressiveStress.plaque(f);
                mediaCompressiveStress(i,j,k)  = R.sumData.compressiveStress.media(f);
                advenCompressiveStress(i,j,k)  = R.sumData.compressiveStress.adventitia(f);
                totalCompressiveStress(i,j,k)  = R.sumData.compressiveStress.total(f);
 
                plaqueEffectiveStrain(i,j,k)   = R.sumData.effectiveStrain.plaque(f);
                mediaEffectiveStrain(i,j,k)    = R.sumData.effectiveStrain.media(f);
                advenEffectiveStrain(i,j,k)    = R.sumData.effectiveStrain.adventitia(f);
                totalEffectiveStrain(i,j,k)    = R.sumData.effectiveStrain.total(f);
 
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
 
% Stenosis
stenVec = repelem(stenosis(:), nK * nR);
 
% Rotation
rotVec  = repmat(rot(:)', nK * nS, 1);
rotVec  = rotVec(:);
 
% Curvature
kapVec  = repmat(kappa(:)', nR * nS, 1);
kapVec  = kapVec(:);
 
nObs_sten = nK * nR;
nObs_rot  = nK * nS;
nObs_kap  = nR * nS;
 
% Build reshaped vectors
rs = @(m) reshapeSet(m, nObs_sten, nS, nObs_kap, nK, nObs_rot, nR);
 
tissueMats = { ...
    struct('plaque', plaqueTensileStress,     'media', mediaTensileStress,     'adven', advenTensileStress), ...
    struct('plaque', plaqueCompressiveStress, 'media', mediaCompressiveStress, 'adven', advenCompressiveStress), ...
    struct('plaque', plaqueTensileStrain,     'media', mediaTensileStrain,     'adven', advenTensileStrain), ...
    struct('plaque', plaqueCompressiveStrain, 'media', mediaCompressiveStrain, 'adven', advenCompressiveStrain), ...
    struct('plaque', plaqueEffectiveStress,   'media', mediaEffectiveStress,   'adven', advenEffectiveStress), ...
    struct('plaque', plaqueEffectiveStrain,   'media', mediaEffectiveStrain,   'adven', advenEffectiveStrain) };
 
measureNames = {'TensileStress','CompressiveStress','TensileStrain','CompressiveStrain','EffectiveStress','EffectiveStrain'};
 
tissues     = {'Plaque', 'Media', 'Adventitia'};
factors     = {'Stenosis', 'Curvature', 'Rotation'};
tissueField = {'plaque','media','adven'};
 
tissueData = cell(1,3);
for ti = 1:3
    s = struct();
    for mi = 1:numel(measureNames)
        s.(measureNames{mi}) = rs(tissueMats{mi}.(tissueField{ti}));
    end
    tissueData{ti} = s;
end
 
factorVecs = {stenVec, kapVec, rotVec};
 
% Compute correlations
plotDefs = { ...
    struct('name','Stress', 'fields',{{'TensileStress','CompressiveStress'}}, 'panelTitles',{{'Tensile','Compressive'}}), ...
    struct('name','Strain', 'fields',{{'TensileStrain','CompressiveStrain'}}, 'panelTitles',{{'Tensile','Compressive'}}), ...
    struct('name','Effective Metrics', 'fields',{{'EffectiveStress','EffectiveStrain'}}, 'panelTitles',{{'Effective Stress','Effective Strain'}}) };
 
rng(42);
for pd = 1:numel(plotDefs)
    def = plotDefs{pd};
    rho_mat  = zeros(3, 3, 2);
    pval_mat = zeros(3, 3, 2);
 
    for ti = 1:3
        for fi = 1:3
            for si = 1:2
                x = factorVecs{fi};
                y = tissueData{ti}.(def.fields{si}){fi};
                [r, p, ~] = spearman_ci(x, y);
                rho_mat(ti,fi,si)  = r;
                pval_mat(ti,fi,si) = p;
            end
        end
    end
 
    figure('Position',[100 100 1500 600], 'Color','white');
    sgtitle(['Spearman \rho Heatmap for ' plottype ' ' def.name], ...
        'FontWeight','bold','FontSize',15);
 
    for si = 1:2
        ax = subplot(1, 2, si);
        rhoSlice = rho_mat(:,:,si);
 
        imagesc(rhoSlice);
        colormap(ax, redblue_cm(256));
        clim([-1 1]);
        colorbar;
 
        set(ax, 'XTick',1:3,'XTickLabel',factors, ...
            'YTick',1:3,'YTickLabel',tissues, ...
            'FontSize',13,'TickDir','out','Box','off');
        title(def.panelTitles{si}, 'FontWeight','bold','FontSize',15);
 
        for ti = 1:3
            for fi = 1:3
                r = rhoSlice(ti,fi);
                textColor = 'white';
                if abs(r) < 0.4, textColor = [0.2 0.2 0.2]; end
                text(fi, ti, sprintf('%.2f', r), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',13,'FontWeight','bold','Color',textColor);
            end
        end
    end
end
 
 
%% Reshape helper
function vecs = reshapeSet(m, nObs_sten, nS, nObs_kap, nK, nObs_rot, nR)
    vSten = reshape(m, nObs_sten, nS);                vSten = vSten(:);
    vKap  = reshape(permute(m,[2 3 1]), nObs_kap, nK); vKap  = vKap(:);
    vRot  = reshape(permute(m,[1 3 2]), nObs_rot, nR); vRot  = vRot(:);
    vecs = {vSten, vKap, vRot};
end
 
%% Red-blue diverging colormap helper
function cm = redblue_cm(n)
    half = floor(n/2);
    r1 = linspace(0.18, 1, half)';
    g1 = linspace(0.45, 1, half)';
    b1 = ones(half, 1);
    r2 = ones(n-half, 1);
    g2 = linspace(1, 0.33, n-half)';
    b2 = linspace(1, 0.10, n-half)';
    cm = [r1 g1 b1; r2 g2 b2];
end 
 
%% Spearman correlation helper
% Returns rho, p-value, and bootstrapped 95% CI
function [rho, pval, ci] = spearman_ci(x, y, nBoot)
    if nargin < 3, nBoot = 2000; end
    keep = (x ~= 0 & ~isnan(x) & y ~= 0 & ~isnan(y));
    x = x(keep); y = y(keep);
    [rho, pval] = corr(x, y, 'Type','Spearman');
    % Bootstrap CI on rho
    rhos = zeros(nBoot, 1);
    n = numel(x);
    for b = 1:nBoot
        idx = randi(n, n, 1);
        [rhos(b), ~] = corr(x(idx), y(idx), 'Type','Spearman');
    end
    ci = quantile(rhos, [0.025 0.975]);
end