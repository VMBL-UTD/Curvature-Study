% Run 96 simulations on HPC using pre-generated meshes

close all
clear
clc

opts                        = struct();
opts.region_length          = 10;
opts.bump_length            = 10;
opts.buffer_length          = 5;
opts.resolution             = 0.18; % 0.15 resolution
opts.lumen_radius           = 1.5;
opts.media_thick            = 0.5;
opts.adventitia_thick       = 0.5;
opts.sleeve_thick           = 0.3;
opts.debug_plots            = 0;
opts.plq_inner_thickness    = 0.3;
inner_plaque_comp           = 'necroticcore';
outer_plaque_comp           = 'necroticcore';
p                           = 0.016;

%% HPC Specific
feb_out_path                = '/groups/hhayenga/udp230000/FEBioStudio';
febio_exe_path              = '/groups/hhayenga/udp230000/FEBioStudio/bin/febio4';
gibbonPath                  = '/groups/hhayenga/udp230000/GIBBON-master';
addpath(genpath(gibbonPath));
codes_path                  = '/groups/hhayenga/udp230000/MATLAB/Tensile Compressive Sims';
addpath(genpath(codes_path));

sav_dir = '/home/udp230000/scratch/HomoNCPlq_Res0.05_NegKappasSims';
mesh_path = '/groups/hhayenga/udp230000/MATLAB/Tensile Compressive Sims/HomoNCPlq_Res0.05_NegKappasSims_Meshes';

%% Sims Loops

% Positive kappas sims for rotation pairs
kappa_list          = [40,80,120];
rotation_list       = [0,45,90,135];
stenosis_list       = [0,0.15,0.3,0.45,0.6,0.75];

for iR = 1:numel(kappa_list)
    k = kappa_list(iR);
    opts.kappa = k;
   for iA = 1:numel(rotation_list)
       r = rotation_list(iA);
       opts.curve_dir = [0, cosd(r), sind(r)];
       for iP = 1:numel(stenosis_list)
           s = stenosis_list(iP);
           opts.percent_stenosis = s;
           
           % build unique file tag
           rTag = sprintf('%+d', k);
           aTag = sprintf('%03d', r);
           sTag = sprintf('%03d', round(s*100));
           feb_fileName = sprintf('k%s_r%s_s%s', rTag, aTag, sTag);
   
           % file locations
           febFile = fullfile(sav_dir, [feb_fileName '.feb']);
           strainCSV = fullfile(sav_dir, [feb_fileName '_strain.csv']);
           stressCSV = fullfile(sav_dir, [feb_fileName '_stress.csv']);
           
           % load pre-generated mesh
           load([mesh_path, '/', feb_fileName, '_data.mat'])

           % write FEBio input
           writeToFEBFile(mesh, inner_plaque_comp, outer_plaque_comp, p, feb_fileName, sav_dir);
           full_path = fullfile(sav_dir, [feb_fileName,'.feb']);
           system(['"',febio_exe_path,'"',' ','"',full_path,'"']);
           
           elem_centroids = getCentroids(mesh.nodes,mesh.elements);
           [axialStrain, radialStrain] = calculateAxialRadialStrain(strainCSV, elem_centroids, mesh, opts);
           [axialStress, radialStress] = calculateAxialRadialStress(stressCSV, elem_centroids, mesh, opts);

           effectiveStress = calculateStress(stressCSV);
           effectiveStrain = calculateStrain(strainCSV);

           axialData = struct('tensileStrain',axialStrain.tensileStrain,'compressiveStrain',axialStrain.compressiveStrain,...
               'tensileStress',axialStress.tensileStress,'compressiveStress',axialStress.compressiveStress);
           effectiveData = struct('effectiveStrain',effectiveStrain,'effectiveStress',effectiveStress);

           dataMat = fullfile(sav_dir, [feb_fileName '_data.mat']);

           % save to .mat files 
           save(dataMat, 'mesh','opts','axialData','effectiveData');

           fprintf('Ran sim for %s\n', feb_fileName);
       end
   end
end

% Negative kappas list for duplicate rotation pairs
kappa_list          = [-40,-80,-120];
rotation_list       = 0;
stenosis_list       = [0,0.15,0.3,0.45,0.6,0.75];

for iR = 1:numel(kappa_list)
    k = kappa_list(iR);
    opts.kappa = k;
   for iA = 1:numel(rotation_list)
       r = rotation_list(iA);
       opts.curve_dir = [0, cosd(r), sind(r)];
       for iP = 1:numel(stenosis_list)
           s = stenosis_list(iP);
           opts.percent_stenosis = s;
           
           % build unique file tag
           rTag = sprintf('%+d', k);
           aTag = sprintf('%03d', r);
           sTag = sprintf('%03d', round(s*100));

           if k < 0
               switch r
                   case 225
                       r = 315;
                   case 270
                       r = 270;
                   case 315
                       r = 225;
                   case 0
                       r = 180;
               end

               feb_fileName = sprintf('k+%d_r%03d_s%s', abs(k), r, sTag);
           end
   
           % file locations
           febFile = fullfile(sav_dir, [feb_fileName '.feb']);
           strainCSV = fullfile(sav_dir, [feb_fileName '_strain.csv']);
           stressCSV = fullfile(sav_dir, [feb_fileName '_stress.csv']);
           
           % load pre-generated mesh
           load([mesh_path, '/', feb_fileName, '_data.mat'])

           % write FEBio input
           writeToFEBFile(mesh, inner_plaque_comp, outer_plaque_comp, p, displacement, febio_exe_path, feb_fileName, sav_dir);
           full_path = fullfile(sav_dir, [feb_fileName,'.feb']);
           system(['"',febio_exe_path,'"',' ','"',full_path,'"']);
           
           elem_centroids = getCentroids(mesh.nodes,mesh.elements);
           axialStrain = calculateAxialStrain(strainCSV, elem_centroids, mesh, opts);
           axialStress = calculateAxialStress(stressCSV, elem_centroids, mesh, opts);

           effectiveStress = calculateEffectiveStress(stressCSV);
           effectiveStrain = calculateEffectiveStrain(strainCSV);

           axialData = struct('tensileStrain',axialStrain.tensileStrain,'compressiveStrain',axialStrain.compressiveStrain,...
               'tensileStress',axialStress.tensileStress,'compressiveStress',axialStress.compressiveStress);
           effectiveData = struct('effectiveStrain',effectiveStrain,'effectiveStress',effectiveStress);

           dataMat = fullfile(sav_dir, [feb_fileName '_data.mat']);

           % save to .mat files 
           save(dataMat, 'mesh','opts','axialData','effectiveData');

           fprintf('Ran sim for %s\n', feb_fileName);
       end
   end
end

% 6 sims for rotation-invariant 0 kappa
kappa_list          = 0;
rotation_list       = 0;
stenosis_list       = [0,0.15,0.3,0.45,0.6,0.75];

for iR = 1:numel(kappa_list)
    k = kappa_list(iR);
    opts.kappa = k;
   for iA = 1:numel(rotation_list)
       r = rotation_list(iA);
       opts.curve_dir = [0, cosd(r), sind(r)];
       for iP = 1:numel(stenosis_list)
           s = stenosis_list(iP);
           opts.percent_stenosis = s;
           
           % build unique file tag
           rTag = sprintf('%+d', k);
           aTag = sprintf('%03d', r);
           sTag = sprintf('%03d', round(s*100));
           feb_fileName = sprintf('k+%d_r%03d_s%s', abs(k), r, sTag);
   
           % file locations
           febFile = fullfile(sav_dir, [feb_fileName '.feb']);
           strainCSV = fullfile(sav_dir, [feb_fileName '_strain.csv']);
           stressCSV = fullfile(sav_dir, [feb_fileName '_stress.csv']);
           
           % load pre-generated mesh
           load([mesh_path, '/', feb_fileName, '_data.mat'])

           % write FEBio input
           writeToFEBFile(mesh, inner_plaque_comp, outer_plaque_comp, p, displacement, febio_exe_path, feb_fileName, sav_dir);
           full_path = fullfile(sav_dir, [feb_fileName,'.feb']);
           system(['"',febio_exe_path,'"',' ','"',full_path,'"']);
           
           elem_centroids = getCentroids(mesh.nodes,mesh.elements);
           axialStrain = calculateAxialStrain(strainCSV, elem_centroids, mesh, opts);
           axialStress = calculateAxialStress(stressCSV, elem_centroids, mesh, opts);

           effectiveStress = calculateEffectiveStress(stressCSV);
           effectiveStrain = calculateEffectiveStrain(strainCSV);

           axialData = struct('tensileStrain',axialStrain.tensileStrain,'compressiveStrain',axialStrain.compressiveStrain,...
               'tensileStress',axialStress.tensileStress,'compressiveStress',axialStress.compressiveStress);
           effectiveData = struct('effectiveStrain',effectiveStrain,'effectiveStress',effectiveStress);

           dataMat = fullfile(sav_dir, [feb_fileName '_data.mat']);

           % save to .mat files 
           save(dataMat, 'mesh','opts','axialData','effectiveData');

           fprintf('Ran sim for %s\n', feb_fileName);
       end
   end
end
