% Generate meshes for running sims on HPC

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

sav_dir = '/groups/hhayenga/udp230000/MATLAB/Tensile Compressive Sims/0.12Meshes';

%% Loops to generate meshes

% Positive kappas sims for rotation pairs
kappa_list          = [40,80,120];
rotation_list       = [0,45,90,135];
stenosis_list       = [0,0.15,0.30,0.45,0.60,0.75];

for iR = 1:numel(kappa_list)
    k = kappa_list(iR);
    opts.kappa = k;
    for iA = 1:numel(rotation_list)
       r = rotation_list(iA);
       opts.curve_dir = [0, cosd(r), sind(r)];
       for iP = 1:numel(stenosis_list)
           s = stenosis_list(iP);
           opts.percent_stenosis = s;
           
           rTag = sprintf('%+d', k);
           aTag = sprintf('%03d', r);
           sTag = sprintf('%03d', round(s*100));
           feb_fileName = sprintf('k%s_r%s_s%s', rTag, aTag, sTag);

           [mesh] = createStenosedArteryMesh(opts);
           dataMat = fullfile(sav_dir, [feb_fileName '_data.mat']);
           
           save(dataMat,'mesh','opts')

           fprintf('Generated mesh for %s\n', feb_fileName);
       end
   end
end

% Negative kappas list for duplicate rotation pairs
kappa_list          = [-40,-80,-120];
rotation_list       = 0;
stenosis_list       = [0,0.15,0.30,0.45,0.60,0.75];

for iR = 1:numel(kappa_list)
    k = kappa_list(iR);
    opts.kappa = k;
    for iA = 1:numel(rotation_list)
       r = rotation_list(iA);
       opts.curve_dir = [0, cosd(r), sind(r)];
       for iP = 1:numel(stenosis_list)
           s = stenosis_list(iP);
           opts.percent_stenosis = s;
           
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

           [mesh] = createStenosedArteryMesh(opts);
           dataMat = fullfile(sav_dir, [feb_fileName '_data.mat']);
           
           save(dataMat,'mesh','opts')

           fprintf('Generated mesh for %s\n', feb_fileName);
       end
   end
end

% 6 meshes for rotation-invariant 0 kappa
kappa_list          = 0;
rotation_list       = 0;
stenosis_list       = [0,0.15,0.30,0.45,0.60,0.75];

for iR = 1:numel(kappa_list)
    k = kappa_list(iR);
    opts.kappa = k;
    for iA = 1:numel(rotation_list)
       r = rotation_list(iA);
       opts.curve_dir = [0, cosd(r), sind(r)];
       for iP = 1:numel(stenosis_list)
           s = stenosis_list(iP);
           opts.percent_stenosis = s;
           
           rTag = sprintf('%+d', k);
           aTag = sprintf('%03d', r);
           sTag = sprintf('%03d', round(s*100));
           feb_fileName = sprintf('k+%d_r%03d_s%s', abs(k), r, sTag);

           [mesh] = createStenosedArteryMesh(opts);
           dataMat = fullfile(sav_dir, [feb_fileName '_data.mat']);
           
           save(dataMat,'mesh','opts')

           fprintf('Generated mesh for %s\n', feb_fileName);
       end
   end
end
