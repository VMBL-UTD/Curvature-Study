function writeToFEBFile(mesh, inner_plaque_comp, outer_plaque_comp, p, feb_fileName, feb_out_path)


media_elements      = mesh.elements(mesh.elementMaterialID == 1, :);
adven_elements      = mesh.elements(mesh.elementMaterialID == 2, :);
sleeve_elements     = mesh.elements(mesh.elementMaterialID == 3, :);
inner_plq_elements  = mesh.elements(mesh.elementMaterialID == 4, :);
outer_plq_elements  = mesh.elements(mesh.elementMaterialID == 5, :);

Fendcaps = mesh.surfs.caps{1}.F;
Flumen = mesh.surfs.lumen.F;

% Write .feb file
febio_spec.ATTR.version                             = '4.0';   % febio_spec version
febio_spec.Module.ATTR.type                         = 'solid'; % Module section

febio_spec.Control.analysis                         = 'STATIC';
febio_spec.Control.time_steps                       = 10;
febio_spec.Control.step_size                        = 0.1;
febio_spec.Control.plot_zero_state                  = 0;
febio_spec.Control.plot_range                       = [0 -1];
febio_spec.Control.plot_stride                      = 1;
febio_spec.Control.adaptor_re_solve                 = 1;

febio_spec.Control.time_stepper.ATTR.type           = 'default';
febio_spec.Control.time_stepper.dtmin               = 0.01;
febio_spec.Control.time_stepper.dtmax               = 0.1;
febio_spec.Control.time_stepper.aggressiveness      = 0;
febio_spec.Control.time_stepper.cutback             = 0.5;
febio_spec.Control.time_stepper.dtforce             = 0;

febio_spec.Control.solver.ATTR.type                 = 'solid';
febio_spec.Control.solver.symmetric_stiffness       = 'symmetric';
febio_spec.Control.solver.equation_scheme           = 'staggered';
febio_spec.Control.solver.equation_order            = 'default';
febio_spec.Control.solver.optimize_bw               = 0;
febio_spec.Control.solver.lstol                     = 0.9;
febio_spec.Control.solver.lsmin                     = 0.01;
febio_spec.Control.solver.lsiter                    = 5;
febio_spec.Control.solver.check_zero_diagonal       = 0;
febio_spec.Control.solver.force_partition           = 0;
febio_spec.Control.solver.zero_diagonal_tol         = 0;
febio_spec.Control.solver.reform_each_time_step     = 1;
febio_spec.Control.solver.reform_augment            = 0;
febio_spec.Control.solver.diverge_reform            = 1;
febio_spec.Control.solver.min_residual              = 1e-20;
febio_spec.Control.solver.max_residual              = 1e20;
febio_spec.Control.solver.dtol                      = 0.001;
febio_spec.Control.solver.etol                      = 0.01;
febio_spec.Control.solver.rtol                      = 0;
febio_spec.Control.solver.rhoi                      = -2;
febio_spec.Control.solver.alpha                     = 1;
febio_spec.Control.solver.beta                      = 0.25;
febio_spec.Control.solver.gamma                     = 0.5;
febio_spec.Control.solver.logSolve                  = 0;
febio_spec.Control.solver.arc_length                = 0;
febio_spec.Control.solver.arc_length_scale          = 0;
febio_spec.Control.solver.qn_method.ATTR.type       = 'BFGS';
febio_spec.Control.solver.qn_method.max_buffer_size = 0;
febio_spec.Control.solver.qn_method.cycle_buffer    = 1;
febio_spec.Control.solver.qn_method.cmax            = 100000;

% GLOBALS
febio_spec.Globals.Constants.T                      = 298;
febio_spec.Globals.Constants.R                      = 8.314e-06;
febio_spec.Globals.Constants.Fc                     = 0;


% MATERIALS
febio_spec.Material.material{1}.ATTR.name           = 'Media Material';
febio_spec.Material.material{1}.ATTR.type           = 'neo-Hookean';
febio_spec.Material.material{1}.ATTR.id             = 1;
febio_spec.Material.material{1}.E                   = 0.3; % Young's Modulus
febio_spec.Material.material{1}.v                   = 0.45; % Poisson's Ratio

febio_spec.Material.material{2}.ATTR.name           = 'Adventitia Material';
febio_spec.Material.material{2}.ATTR.type           = 'neo-Hookean';
febio_spec.Material.material{2}.ATTR.id             = 2;
febio_spec.Material.material{2}.E                   = 0.8;
febio_spec.Material.material{2}.v                   = 0.45;

febio_spec.Material.material{3}.ATTR.name           = 'Sleeve Material';
febio_spec.Material.material{3}.ATTR.type           = 'neo-Hookean';
febio_spec.Material.material{3}.ATTR.id             = 3;
febio_spec.Material.material{3}.E                   = 0.4;
febio_spec.Material.material{3}.v                   = 0.48;

if mesh.stenosis == 1
    febio_spec.Material.material{4}.ATTR.name           = 'Inner Plaque Material';
    febio_spec.Material.material{4}.ATTR.type           = 'neo-Hookean';
    febio_spec.Material.material{4}.ATTR.id             = 4;
    if strcmp(inner_plaque_comp, 'fibrotic')
        febio_spec.Material.material{4}.E                   = 0.6;
        febio_spec.Material.material{4}.v                   = 0.48;
    elseif strcmp(inner_plaque_comp, 'fibrofatty')
        febio_spec.Material.material{4}.E                   = 0.5;
        febio_spec.Material.material{4}.v                   = 0.48;
    elseif strcmp(inner_plaque_comp, 'calcium')
        febio_spec.Material.material{4}.E                   = 10;
        febio_spec.Material.material{4}.v                   = 0.48;
    elseif strcmp(inner_plaque_comp, 'necroticcore')
        febio_spec.Material.material{4}.E                   = 0.02;
        febio_spec.Material.material{4}.v                   = 0.48;
    elseif strcmp(inner_plaque_comp,'media')
        febio_spec.Material.material{4}.E                   = 0.3;
        febio_spec.Material.material{4}.v                   = 0.45;
    end
    
    febio_spec.Material.material{5}.ATTR.name           = 'Outer Plaque Material';
    febio_spec.Material.material{5}.ATTR.type           = 'neo-Hookean';
    febio_spec.Material.material{5}.ATTR.id             = 5;
    if strcmp(outer_plaque_comp, 'fibrotic')  
        febio_spec.Material.material{5}.E                   = 0.6;
        febio_spec.Material.material{5}.v                   = 0.48;
    elseif strcmp(outer_plaque_comp, 'fibrofatty')
        febio_spec.Material.material{5}.E                   = 0.5;
        febio_spec.Material.material{5}.v                   = 0.48;
    elseif strcmp(outer_plaque_comp, 'calcium')
        febio_spec.Material.material{5}.E                   = 10;
        febio_spec.Material.material{5}.v                   = 0.48;
    elseif strcmp(outer_plaque_comp, 'necroticcore')
        febio_spec.Material.material{5}.E                   = 0.02;
        febio_spec.Material.material{5}.v                   = 0.48;
    elseif strcmp(outer_plaque_comp,'media')
        febio_spec.Material.material{5}.E                   = 0.3;
        febio_spec.Material.material{5}.v                   = 0.45;
    end
end

%% MESH
% 
media_ind      = find(mesh.elementMaterialID == 1);
adven_ind      = find(mesh.elementMaterialID == 2);
sleeve_ind      = find(mesh.elementMaterialID == 3);
inner_plq_ind   = find(mesh.elementMaterialID == 4);
outer_plq_ind   = find(mesh.elementMaterialID == 5);


febio_spec.Mesh.Nodes{1}.ATTR.name                  = 'allNodes'; %The node set name
febio_spec.Mesh.Nodes{1}.node.ATTR.id               = (1:1:size(mesh.nodes,1))'; %The node id's
febio_spec.Mesh.Nodes{1}.node.VAL                   = mesh.nodes; %The nodal coordinates

% Media
febio_spec.Mesh.Elements{1}.ATTR.type               = 'tet10';
febio_spec.Mesh.Elements{1}.ATTR.mat                = 1;
febio_spec.Mesh.Elements{1}.ATTR.name               = 'Media Elements';
febio_spec.Mesh.Elements{1}.elem.ATTR.id            = media_ind; %(1:1:length(media_ind))';
febio_spec.Mesh.Elements{1}.elem.VAL                = media_elements;

febio_spec.MeshDomains.SolidDomain{1}.ATTR.name     = 'Media Elements';
febio_spec.MeshDomains.SolidDomain{1}.ATTR.mat      = 'Media Material';

% Adventitia
febio_spec.Mesh.Elements{2}.ATTR.type               = 'tet10';
febio_spec.Mesh.Elements{2}.ATTR.mat                = 2;
febio_spec.Mesh.Elements{2}.ATTR.name               = 'Adventitia Elements';
febio_spec.Mesh.Elements{2}.elem.ATTR.id            = adven_ind;  %(1:1:length(adven_ind))' + length(media_ind);
febio_spec.Mesh.Elements{2}.elem.VAL                = adven_elements;

febio_spec.MeshDomains.SolidDomain{2}.ATTR.name     = 'Adventitia Elements';
febio_spec.MeshDomains.SolidDomain{2}.ATTR.mat      = 'Adventitia Material';

% Sleeve
febio_spec.Mesh.Elements{3}.ATTR.type               = 'tet10';
febio_spec.Mesh.Elements{3}.ATTR.mat                = 3;
febio_spec.Mesh.Elements{3}.ATTR.name               = 'Sleeve Elements';
febio_spec.Mesh.Elements{3}.elem.ATTR.id            = sleeve_ind; %(1:1:length(sleeve_ind))' + length(adven_ind) + length(media_ind); % sleeve_ind
febio_spec.Mesh.Elements{3}.elem.VAL                = sleeve_elements;

febio_spec.MeshDomains.SolidDomain{3}.ATTR.name     = 'Sleeve Elements';
febio_spec.MeshDomains.SolidDomain{3}.ATTR.mat      = 'Sleeve Material';

if mesh.stenosis == 1
    % Inner Plaque
    febio_spec.Mesh.Elements{4}.ATTR.type               = 'tet10';
    febio_spec.Mesh.Elements{4}.ATTR.mat                = 4;
    febio_spec.Mesh.Elements{4}.ATTR.name               = 'Inner Plaque Elements';
    febio_spec.Mesh.Elements{4}.elem.ATTR.id            = inner_plq_ind; %(1:1:length(inner_plq_ind))' + length(sleeve_ind) + length(adven_ind) + length(media_ind);
    febio_spec.Mesh.Elements{4}.elem.VAL                = inner_plq_elements;
    
    febio_spec.MeshDomains.SolidDomain{4}.ATTR.name     = 'Inner Plaque Elements';
    febio_spec.MeshDomains.SolidDomain{4}.ATTR.mat      = 'Inner Plaque Material';
    
    % Outer Plaque
    febio_spec.Mesh.Elements{5}.ATTR.type               = 'tet10';
    febio_spec.Mesh.Elements{5}.ATTR.mat                = 5;
    febio_spec.Mesh.Elements{5}.ATTR.name               = 'Outer Plaque Elements';
    febio_spec.Mesh.Elements{5}.elem.ATTR.id            = outer_plq_ind; %(1:1:length(outer_plq_ind))' + length(inner_plq_ind) + length(sleeve_ind) + length(adven_ind) + length(media_ind);
    febio_spec.Mesh.Elements{5}.elem.VAL                = outer_plq_elements;
    
    febio_spec.MeshDomains.SolidDomain{5}.ATTR.name     = 'Outer Plaque Elements';
    febio_spec.MeshDomains.SolidDomain{5}.ATTR.mat      = 'Outer Plaque Material';
end

%% DEFINE SURFACES
febio_spec.Mesh.Surface{1}.ATTR.name                = 'Fixed End Surfaces';
febio_spec.Mesh.Surface{1}.tri6.VAL                 = mesh.facesBoundary(Fendcaps,[1,3,5,2,4,6]);
febio_spec.Mesh.Surface{1}.tri6.ATTR.id             = Fendcaps;

febio_spec.Mesh.Surface{2}.ATTR.name                = 'Lumen Surface';
febio_spec.Mesh.Surface{2}.tri6.VAL                 = mesh.facesBoundary(Flumen,[1,3,5,2,4,6]);
febio_spec.Mesh.Surface{2}.tri6.ATTR.id             = Flumen;

%% BOUNDARY
febio_spec.Boundary.bc{1}.ATTR.name                     = 'Fixed Ends';
febio_spec.Boundary.bc{1}.ATTR.node_set                 = strcat('@surface:','Fixed End Surfaces');
febio_spec.Boundary.bc{1}.ATTR.type                     = 'zero displacement';
febio_spec.Boundary.bc{1}.x_dof                         = 1;
febio_spec.Boundary.bc{1}.y_dof                         = 1;
febio_spec.Boundary.bc{1}.z_dof                         = 1;

%% LOADS
febio_spec.Loads.surface_load.ATTR.name                 = 'Lumen Pressure';
febio_spec.Loads.surface_load.ATTR.surface              = 'Lumen Surface';
febio_spec.Loads.surface_load.ATTR.type                 = 'pressure';
febio_spec.Loads.surface_load.pressure.ATTR.lc          = '1'; % load controller id
febio_spec.Loads.surface_load.pressure.VAL              = p;
febio_spec.Loads.surface_load.symmetric_stiffness       = 1;
febio_spec.Loads.surface_load.linear                    = 0;
febio_spec.Loads.surface_load.shell_bottom              = 0;

febio_spec.LoadData.load_controller{1}.ATTR.name           = 'Pressurize';
febio_spec.LoadData.load_controller{1}.ATTR.id             = 1;
febio_spec.LoadData.load_controller{1}.ATTR.type           = 'loadcurve';
febio_spec.LoadData.load_controller{1}.interpolate         = 'LINEAR';
febio_spec.LoadData.load_controller{1}.points.point.VAL    = [0 0; 1 1];

%% LOG FILES

febio_spec.Output.logfile.element_data{1}.ATTR.file        = strcat(feb_fileName, '_strain.csv');
febio_spec.Output.logfile.element_data{1}.ATTR.data        = 'Ex;Ey;Ez;Exy;Eyz;Exz';
febio_spec.Output.logfile.element_data{1}.ATTR.delim       = ',';

febio_spec.Output.logfile.element_data{2}.ATTR.file        = strcat(feb_fileName, '_stress.csv');
febio_spec.Output.logfile.element_data{2}.ATTR.data        = 'sx;sy;sz;sxy;syz;sxz';
febio_spec.Output.logfile.element_data{2}.ATTR.delim       = ',';


%% Generate .feb file;
full_path = fullfile(feb_out_path, [feb_fileName,'.feb']);
disp(full_path)

% Generate .feb file;
febioStruct2xml(febio_spec,full_path);

end