% Post processing script
% Calculates mean or max effective, tensile & compressive
% stress and strain for each simulation
% for elements only in region of interest (ROI)

%% Resolution 0.12 mm Post Processing
clc
clear

feb_out_path_012 = "F:\NCUR Simulations\8.4_0.12Res";
save_path_012    = "F:\NCUR Simulations\8.4_0.12Res";
calcType = 'Max';

results012 = batchProcessFolder(feb_out_path_012, calcType);

dataMat012 = fullfile(save_path_012, [lower(calcType) 'CalcResults_res012.mat']);
save(dataMat012, '-struct', 'results012');

fprintf('Saved %s\n', dataMat012);

%% Resolution 0.10 mm Post Processing
clc
clear

feb_out_path_010 = "F:\NCUR Simulations\4.2_0.1Res_NCPlaque";
save_path_010    = "F:\NCUR Simulations\4.2_0.1Res_NCPlaque";
calcType = 'Max';

results010 = batchProcessFolder(feb_out_path_010, calcType);

dataMat010 = fullfile(save_path_010, [lower(calcType) 'CalcResults_res010.mat']);
save(dataMat010, '-struct', 'results010');

fprintf('Saved %s\n', dataMat010);

%% Helper functions

function results = batchProcessFolder(feb_out_path, calcType)
% Loops over every *_data.mat file in feb_out_path and aggregates
% per-sim results into arrays using processSimEffective.

    files = dir(fullfile(feb_out_path, '*_data.mat'));
    nRuns = numel(files);

    kappa       = zeros(nRuns,1);
    rotation    = zeros(nRuns,1);
    stenosis    = zeros(nRuns,1);
    resolution  = zeros(nRuns,1);
    edge_length = zeros(nRuns,1);

    template = struct('total',      zeros(nRuns,1), ...
                       'plaque',     zeros(nRuns,1), ...
                       'media',      zeros(nRuns,1), ...
                       'adventitia', zeros(nRuns,1));

    sumData = struct( ...
        'effectiveStrain', template, ...
        'effectiveStress', template, ...
        'tensileStrain', template, ...
        'tensileStress', template, ...
        'compressiveStrain', template, ...
        'compressiveStress', template);

    tissueFields = {'total', 'plaque', 'media', 'adventitia'};
    quantFields  = {'effectiveStrain', 'effectiveStress', ...
                    'tensileStrain', 'tensileStress',...
                    'compressiveStrain', 'compressiveStress'};

    for i = 1:nRuns
        fname = files(i).name;

        [k_i, r_i, s_i, res_i, el_i, sum_i] = processSim(feb_out_path, fname, calcType);

        kappa(i)       = k_i;
        rotation(i)    = r_i;
        stenosis(i)    = s_i;
        resolution(i)  = res_i;
        edge_length(i) = el_i;

        for qi = 1:length(quantFields)
            for ti = 1:length(tissueFields)
                tf = tissueFields{ti};
                qf = quantFields{qi};
                if isfield(sum_i.(qf), tf)
                    sumData.(qf).(tf)(i) = sum_i.(qf).(tf);
                else
                    sumData.(qf).(tf)(i) = NaN;
                end
            end
        end
    end

    results.kappa            = kappa;
    results.rotation         = rotation;
    results.stenosis         = stenosis;
    results.resolution       = resolution;
    results.edge_length      = edge_length;
    results.sumData          = sumData;
end

function [kappa, rotation, stenosis, resolution, edge_length, ...
          sumData] = processSim(feb_out_path, fname, calcType)

    MAT_Media = 1;
    MAT_Adven = 2;

    % Parse kappa, rotation, stenosis from filename
    tokens = regexp(fname, 'k([+-]\d+)_r(\d+)_s(\d+)', 'tokens');
    tokens = tokens{1};
    k_val = str2double(tokens{1});
    r_val = str2double(tokens{2});
    s_val = str2double(tokens{3}) / 100;

    L = load(fullfile(feb_out_path, fname), 'mesh', 'opts', 'effectiveData', 'axialData');
    effStrain = L.effectiveData.effectiveStrain;
    effStress = L.effectiveData.effectiveStress;
    tensileStress = L.axialData.tensileStress;
    tensileStrain = L.axialData.tensileStrain;
    compressiveStress = L.axialData.compressiveStress;
    compressiveStrain = L.axialData.compressiveStrain;

    mesh = L.mesh;
    opts = L.opts;
    opts.percent_stenosis = s_val;

    if r_val == 180
        opts.kappa     = -abs(k_val);
        opts.curve_dir = [0, 1, 0];
    elseif r_val == 315
        opts.kappa     = -abs(k_val);
        opts.curve_dir = [0, cosd(225), sind(225)];
    elseif r_val == 270
        opts.kappa     = -abs(k_val);
        opts.curve_dir = [0, cosd(270), sind(270)];
    elseif r_val == 225
        opts.kappa     = -abs(k_val);
        opts.curve_dir = [0, cosd(315), sind(315)];
    else
        opts.kappa     = abs(k_val);
        opts.curve_dir = [0, cosd(r_val), sind(r_val)];
    end

    kappa          = opts.kappa;
    stenosis       = opts.percent_stenosis;
    rotation       = atan2d(opts.curve_dir(3), opts.curve_dir(2));
    resolution     = opts.resolution;
    edge_length    = mean(patchEdgeLengths(mesh.elements, mesh.nodes));

    % Masks
    maskROI = calculateROI(mesh, opts);
    emID    = mesh.elementMaterialID;

    % Total (ROI)
    if strcmp(calcType,"Max")
        sumData.effectiveStrain.total = maxNonzero(effStrain, maskROI);
        sumData.effectiveStress.total = maxNonzero(effStress, maskROI);
        sumData.tensileStrain.total = maxNonzero(tensileStrain, maskROI);
        sumData.tensileStress.total = maxNonzero(tensileStress, maskROI);
        sumData.compressiveStrain.total = maxNonzero(compressiveStrain, maskROI);
        sumData.compressiveStress.total = maxNonzero(compressiveStress, maskROI);
    else
        sumData.effectiveStrain.total = meanNonzero(effStrain, maskROI);
        sumData.effectiveStress.total = meanNonzero(effStress, maskROI);
        sumData.tensileStrain.total = meanNonzero(tensileStrain, maskROI);
        sumData.tensileStress.total = meanNonzero(tensileStress, maskROI);
        sumData.compressiveStrain.total = meanNonzero(compressiveStrain, maskROI);
        sumData.compressiveStress.total = meanNonzero(compressiveStress, maskROI);
    end

    % Tissue-specific (ROI)
    materials = {
        [],        "plaque";
        MAT_Media, "media";
        MAT_Adven, "adventitia"
    };

    for m = 1:size(materials,1)
        label = materials{m,2};

        if strcmp(label, 'plaque')
            maskM = maskROI & (emID == 4 | emID == 5);
        else
            matID = materials{m,1};
            maskM = maskROI & (emID == matID);
        end

        if any(maskM)
            if strcmp(calcType, 'Max')
                sumData.effectiveStrain.(label) = maxNonzero(effStrain, maskM);
                sumData.effectiveStress.(label) = maxNonzero(effStress, maskM);
                sumData.tensileStrain.(label) = maxNonzero(tensileStrain, maskM);
                sumData.tensileStress.(label) = maxNonzero(tensileStress, maskM);
                sumData.compressiveStrain.(label) = maxNonzero(compressiveStrain, maskM);
                sumData.compressiveStress.(label) = maxNonzero(compressiveStress, maskM);
            else
                sumData.effectiveStrain.(label) = NaN;
                sumData.effectiveStress.(label) = NaN;
                sumData.tensileStrain.(label) = NaN;
                sumData.tensileStress.(label) = NaN;
                sumData.compressiveStrain.(label) = NaN;
                sumData.compressiveStress.(label) = NaN;
            end
        else
            if strcmp(calcType, 'Max')
                sumData.effectiveStrain.(label) = meanNonzero(effStrain, maskM);
                sumData.effectiveStress.(label) = meanNonzero(effStress, maskM);
                sumData.tensileStrain.(label) = meanNonzero(tensileStrain, maskM);
                sumData.tensileStress.(label) = meanNonzero(tensileStress, maskM);
                sumData.compressiveStrain.(label) = meanNonzero(compressiveStrain, maskM);
                sumData.compressiveStress.(label) = meanNonzero(compressiveStress, maskM);
            else
                sumData.effectiveStrain.(label) = NaN;
                sumData.effectiveStress.(label) = NaN;
                sumData.tensileStrain.(label) = NaN;
                sumData.tensileStress.(label) = NaN;
                sumData.compressiveStrain.(label) = NaN;
                sumData.compressiveStress.(label) = NaN;
            end
        end
    end

    fprintf('Processed %s | kappa=%.0f rotation=%.0f stenosis=%.0f%%\n', ...
        fname, kappa, rotation, stenosis*100);
end

function m = maxNonzero(X, mask)
    vals = X(mask, end);
    vals = vals(vals ~= 0);
    if isempty(vals)
        m = NaN;
    else
        if mean(vals) < 0
            m = min(vals);  % most negative = peak compression
        else
            m = max(vals);  % most positive = peak tension
        end
    end
end

function m = meanNonzero(X, mask)
    vals = X(mask, end);
    vals = vals(vals ~= 0);
    if isempty(vals)
        m = NaN;
    else
        m = mean(vals);
    end
end
