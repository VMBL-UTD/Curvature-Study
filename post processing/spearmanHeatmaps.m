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

% Fill matrices
for i = 1:nK
    for j = 1:nR
        for k = 1:nS
            idx = (abs(R.kappa) == kappa(i)) & ...
                  (wrapTo360(R.rotation) == rot(j)) & ...
                  (R.stenosis*100 == stenosis(k));
            if any(idx)
                f = find(idx,1,'first');
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

% Looping for stress & strain
for i=1:2
    switch i
        case 1
            value = 'Stress';
            plaqueTensile = plaqueTensileStress;
            plaqueCompressive = plaqueCompressiveStress;
            mediaTensile = mediaTensileStress;
            mediaCompressive = mediaCompressiveStress;
            advenTensile = advenTensileStress;
            advenCompressive = advenCompressiveStress;
        otherwise
            value = 'Strain';
            plaqueTensile = plaqueTensileStrain;
            plaqueCompressive = plaqueCompressiveStrain;
            mediaTensile = mediaTensileStrain;
            mediaCompressive = mediaCompressiveStrain;
            advenTensile = advenTensileStrain;
            advenCompressive = advenCompressiveStrain;
    end
    % Tensile
    pT_sten = reshape(plaqueTensile, nObs_sten, nS); pT_sten = pT_sten(:);
    mT_sten = reshape(mediaTensile, nObs_sten, nS); mT_sten = mT_sten(:);
    aT_sten = reshape(advenTensile, nObs_sten, nS); aT_sten = aT_sten(:);

    pT_rot  = reshape(permute(plaqueTensile, [1 3 2]), nObs_rot, nR); pT_rot = pT_rot(:);
    mT_rot  = reshape(permute(mediaTensile, [1 3 2]), nObs_rot, nR); mT_rot = mT_rot(:);
    aT_rot  = reshape(permute(advenTensile, [1 3 2]), nObs_rot, nR); aT_rot = aT_rot(:);

    pT_kap  = reshape(permute(plaqueTensile, [2 3 1]), nObs_kap, nK); pT_kap = pT_kap(:);
    mT_kap  = reshape(permute(mediaTensile, [2 3 1]), nObs_kap, nK); mT_kap = mT_kap(:);
    aT_kap  = reshape(permute(advenTensile, [2 3 1]), nObs_kap, nK); aT_kap = aT_kap(:);

    % Compressive
    pC_sten = reshape(plaqueCompressive, nObs_sten, nS); pC_sten = pC_sten(:);
    mC_sten = reshape(mediaCompressive, nObs_sten, nS); mC_sten = mC_sten(:);
    aC_sten = reshape(advenCompressive, nObs_sten, nS); aC_sten = aC_sten(:);

    pC_rot  = reshape(permute(plaqueCompressive, [1 3 2]), nObs_rot, nR); pC_rot = pC_rot(:);
    mC_rot  = reshape(permute(mediaCompressive, [1 3 2]), nObs_rot, nR); mC_rot = mC_rot(:);
    aC_rot  = reshape(permute(advenCompressive, [1 3 2]), nObs_rot, nR); aC_rot = aC_rot(:);

    pC_kap  = reshape(permute(plaqueCompressive, [2 3 1]), nObs_kap, nK); pC_kap = pC_kap(:);
    mC_kap  = reshape(permute(mediaCompressive, [2 3 1]), nObs_kap, nK); mC_kap = mC_kap(:);
    aC_kap  = reshape(permute(advenCompressive, [2 3 1]), nObs_kap, nK); aC_kap = aC_kap(:);

    % Compute correlations
    tissues  = {'Plaque', 'Media', 'Adventitia'};
    factors  = {'Stenosis', 'Curvature', 'Rotation'};
    stypes   = {'Tensile', 'Compressive'};

    rho_mat   = zeros(3, 3, 2);
    pval_mat  = zeros(3, 3, 2);
    ci_lo_mat = zeros(3, 3, 2);
    ci_hi_mat = zeros(3, 3, 2);

    % fi=1 → Stenosis, fi=2 → Curvature, fi=3 → Rotation
    stressVecs{1}{1}{1} = pT_sten; stressVecs{1}{1}{2} = pC_sten;
    stressVecs{1}{2}{1} = pT_kap;  stressVecs{1}{2}{2} = pC_kap;
    stressVecs{1}{3}{1} = pT_rot;  stressVecs{1}{3}{2} = pC_rot;

    stressVecs{2}{1}{1} = mT_sten; stressVecs{2}{1}{2} = mC_sten;
    stressVecs{2}{2}{1} = mT_kap;  stressVecs{2}{2}{2} = mC_kap;
    stressVecs{2}{3}{1} = mT_rot;  stressVecs{2}{3}{2} = mC_rot;

    stressVecs{3}{1}{1} = aT_sten; stressVecs{3}{1}{2} = aC_sten;
    stressVecs{3}{2}{1} = aT_kap;  stressVecs{3}{2}{2} = aC_kap;
    stressVecs{3}{3}{1} = aT_rot;  stressVecs{3}{3}{2} = aC_rot;

    factorVecs = {stenVec, kapVec, rotVec};

    rng(42);
    for ti = 1:3
        for fi = 1:3
            for si = 1:2
                x = factorVecs{fi};
                y = stressVecs{ti}{fi}{si};
                [r, p, ci] = spearman_ci(x, y);
                rho_mat(ti,fi,si)   = r;
                pval_mat(ti,fi,si)  = p;
                ci_lo_mat(ti,fi,si) = ci(1);
                ci_hi_mat(ti,fi,si) = ci(2);
            end
        end
    end

    % Summary Heatmap
    % Cell value = rho, color intensity = |rho|, red/blue = sign

    figure('Position',[100 100 1500 600], 'Color','white', ...
        'Name','Sensitivity heatmap — all tissues');
    sgtitle(['Spearman \rho Heatmap for ' plottype ' ' value], ...
        'FontWeight','bold','FontSize',15);

    for si = 1:2
        ax = subplot(1, 2, si);
        rhoSlice = rho_mat(:,:,si); % 3 tissues × 3 factors

        imagesc(rhoSlice);
        colormap(ax, redblue_cm(256));
        clim([-1 1]);
        colorbar;

        set(ax, 'XTick',1:3,'XTickLabel',factors, ...
            'YTick',1:3,'YTickLabel',tissues, ...
            'FontSize',13,'TickDir','out','Box','off');
        title(stypes{si},'FontWeight','bold','FontSize',15);

        % Annotate cells
        for ti = 1:3
            for fi = 1:3
                r = rhoSlice(ti,fi);
                textColor = 'white';
                if abs(r) < 0.4, textColor = [0.2 0.2 0.2]; end
                text(fi, ti, sprintf('%.2f%s', r), ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'FontSize',13,'FontWeight','bold','Color',textColor);
            end
        end
    end
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