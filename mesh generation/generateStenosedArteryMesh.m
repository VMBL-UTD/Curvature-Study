%{

Function to generate a volumetric tetrahedral mesh of a stenosed artery. Returns a struct object with the following fields:
    - 'nodes' Array containing all mesh nodes
    - 'facesBoundary' Array matrix containing the boundary faces of the mesh
    - 'faces' Array containing the faces of all elements
    - 'elements' Array containing the element connectivity
    - 'elementMaterialID' Array with the material id of each element (1: adventitia, 2: media, 7: intima, 8: sleeve)
    - 'faceMaterialID' Array corresponding to the 'faces' array that defines the material of each face (for plotting)
    - 'surfs' Struct containing surface definitions for the lumen, exterior, and end cap surfaces
    - 'ec' Array containing the centroids of each element
    - 'element_volume' Array containing the volume of each element
    - 'inletNodes' Nodal indices on the inlet edge
    - 'outletNodes' Nodal indices on the outlet edge
Plaque is not automatically seeded in the mesh, must be added after mesh is generated.


Example use case for an artery with a 50% obstruction
opts = struct();
opts.region_length = 10;
opts.bump_length = 10;
opts.buffer_length = 5;
opts.resolution = 0.26;
opts.lumen_radius = 1.5;
opts.media_thick = 0.5;
opts.adventitia_thick = 0.5;
opts.sleeve_thick = 0.3;
opts.percent_stenosis = 0.5;
opts.debug_plots = 0;

opts.plq_inner_thickness = 0.3;
opts.kappa = 80;
opts.curve_dir = [0,1,1];

[mesh] = createStenosedArtery_wplaque2(opts);

%}

function [mesh] = generateStenosedArteryMesh(opts)

if isfield(opts,'region_length'); region_length = opts.region_length; else; error("Field 'region_length' missing from input struct"); end
if isfield(opts,'bump_length'); bump_length = opts.bump_length; else; error("Field 'bump_length' missing from input struct"); end
if isfield(opts,'buffer_length'); buffer_length = opts.buffer_length; else; error("Field 'buffer_length' missing from input struct"); end
if isfield(opts,'resolution'); resolution = opts.resolution; else; error("Field 'resolution' missing from input struct"); end
if isfield(opts,'lumen_radius'); lumen_radius = opts.lumen_radius; else; error("Field 'lumen_radius' missing from input struct"); end
if isfield(opts,'media_thick'); media_thick = opts.media_thick; else; error("Field 'media_thick' missing from input struct"); end
if isfield(opts,'adventitia_thick'); adventitia_thick = opts.adventitia_thick; else; error("Field 'adventitia_thick' missing from input struct"); end
if isfield(opts,'sleeve_thick'); sleeve_thick = opts.sleeve_thick; else; error("Field 'sleeve_thick' missing from input struct"); end
if isfield(opts,'percent_stenosis'); percent_stenosis = opts.percent_stenosis; else; error("Field 'percent_stenosis' missing from input struct"); end
if isfield(opts,'kappa'); kappa = opts.kappa; else; error("Field 'kappa' missing from input struct"); end
if isfield(opts,'debug_plots'); debug_plots = opts.debug_plots; else; debug_plots = 0; end
if isfield(opts,'plq_inner_thickness'); plq_inner_thickness = opts.plq_inner_thickness; else; error("Field 'plq_inner_thickness' missing from input struct"); end
if isfield(opts,'curve_dir'); curve_dir = opts.curve_dir; else; error("Field 'curve_dir' missing from input struct"); end

curve_dir = curve_dir./norm(curve_dir);


%%
% Define cylinder mesh parameters
num_radial = floor(2*pi*lumen_radius/resolution);
num_height = floor(region_length/(resolution*0.866));

optstr.pre.max_hole_area = 100;
optstr.pre.max_hole_edges = 0;
optstr.pointSpacing = resolution;

optionStruct.outputType='label';



%% Create lumen mesh
[F1,V1,v1top,v1bot] = createStenosedLumen(opts);
F1 = fliplr(F1);
p1 = V1(v1top,:);
p2 = V1(v1bot,:);



%% Create media->adventitia boundary
clen = region_length+(2*buffer_length);
[F2,Vt] = patchcylinder(lumen_radius+media_thick,num_radial,clen,num_height,'tri');
V2 = fliplr(Vt) + [clen/2 0 0];
V2 = applyCurvature(V2,kappa,curve_dir);
[F2,V2] = ggremesh(F2,V2,optstr);

% Get mesh edge boundaries
[E,~,~] = patchBoundary(F2);
optionStruct.outputType='label';
g1=tesgroup(E,optionStruct);
ve1 = edgeListToCurve(E(g1==1,:));
ve2 = edgeListToCurve(E(g1==2,:));

d1 = sqrt(sum(mean(V2(ve1,:),1).^2,2));
d2 = sqrt(sum(mean(V2(ve2,:),1).^2,2));

if d1>d2
    v2top = ve2;
    v2bot = ve1;
else
    v2top = ve1;
    v2bot = ve2;
end

pd1 = zeros(length(v2top),1);
for i=1:length(v2top)
    pd1(i,:) = pointToPlaneDistance(mean(p1,1),p1(2,:),p1(3,:),V2(v2top(i),:));
end

pd2 = zeros(length(v2bot),1);
for i=1:length(v2bot)
    pd2(i,:) = pointToPlaneDistance(mean(p2,1),p2(2,:),p2(3,:),V2(v2bot(i),:));
end

V2(v2top,:) = V2(v2top,:)-[1,0,0].*pd1;
V2(v2bot,:) = V2(v2bot,:)-cross(p2(1,:)-mean(p2,1),p2(2,:)-mean(p2,1)).*pd2;



%% Create adventitia->sleeve boundary
clen = region_length+(2*buffer_length);
[F3,Vt] = patchcylinder(lumen_radius+media_thick+adventitia_thick,num_radial,clen,num_height,'tri');
V3 = fliplr(Vt) + [clen/2 0 0];
V3 = applyCurvature(V3,kappa,curve_dir);
[F3,V3] = ggremesh(F3,V3,optstr);

% Get mesh edge boundaries
[E,~,~] = patchBoundary(F3);
optionStruct.outputType='label';
g1=tesgroup(E,optionStruct);
ve1 = edgeListToCurve(E(g1==1,:));
ve2 = edgeListToCurve(E(g1==2,:));

d1 = sqrt(sum(mean(V3(ve1,:),1).^2,2));
d2 = sqrt(sum(mean(V3(ve2,:),1).^2,2));

if d1>d2
    v3top = ve2;
    v3bot = ve1;
else
    v3top = ve1;
    v3bot = ve2;
end

pd1 = zeros(length(v3top),1);
for i=1:length(v3top)
    pd1(i,:) = pointToPlaneDistance(mean(p1,1),p1(2,:),p1(3,:),V3(v3top(i),:));
end

pd2 = zeros(length(v3bot),1);
for i=1:length(v3bot)
    pd2(i,:) = pointToPlaneDistance(mean(p2,1),p2(2,:),p2(3,:),V3(v3bot(i),:));
end

V3(v3top,:) = V3(v3top,:)-[1,0,0].*pd1;
V3(v3bot,:) = V3(v3bot,:)-cross(p2(1,:)-mean(p2,1),p2(2,:)-mean(p2,1)).*pd2;



%% Create outer mesh
clen = region_length+(2*buffer_length);
[F4,Vt] = patchcylinder(lumen_radius+media_thick+adventitia_thick+sleeve_thick,num_radial,clen,num_height,'tri');
V4 = fliplr(Vt) + [clen/2 0 0];
V4 = applyCurvature(V4,kappa,curve_dir);
[F4,V4] = ggremesh(F4,V4,optstr);

% Get mesh edge boundaries
[E,~,~] = patchBoundary(F4);
optionStruct.outputType='label';
g1=tesgroup(E,optionStruct);
ve1 = edgeListToCurve(E(g1==1,:));
ve2 = edgeListToCurve(E(g1==2,:));

d1 = sqrt(sum(mean(V4(ve1,:),1).^2,2));
d2 = sqrt(sum(mean(V4(ve2,:),1).^2,2));

if d1>d2
    v4top = ve2;
    v4bot = ve1;
else
    v4top = ve1;
    v4bot = ve2;
end

pd1 = zeros(length(v4top),1);
for i=1:length(v4top)
    pd1(i,:) = pointToPlaneDistance(mean(p1,1),p1(2,:),p1(3,:),V4(v4top(i),:));
end

pd2 = zeros(length(v4bot),1);
for i=1:length(v4bot)
    pd2(i,:) = pointToPlaneDistance(mean(p2,1),p2(2,:),p2(3,:),V4(v4bot(i),:));
end

V4(v4top,:) = V4(v4top,:)-[1,0,0].*pd1;
V4(v4bot,:) = V4(v4bot,:)-cross(p2(1,:)-mean(p2,1),p2(2,:)-mean(p2,1)).*pd2;



if debug_plots
    figure
    hold on

    gpatch(F1,V1,'r');
    plot3V(V1(v1top,:),'r*')
    plot3V(V1(v1bot,:),'b*')

    gpatch(F2,V2,'b');
    plot3V(V2(v2top,:),'r*')
    plot3V(V2(v2bot,:),'b*')

    gpatch(F3,V3,'g');
    plot3V(V3(v3top,:),'r*')
    plot3V(V3(v3bot,:),'b*')

    gpatch(F4,V4,'m');
    plot3V(V4(v4top,:),'r*')
    plot3V(V4(v4bot,:),'b*')

    legend({'Lumen','','','Med-Adv','','','Adv-Slv','','','Outer'})
    axisGeom
end


%%
clc
[Ft1,Vt1] = regionTriMesh3D({V1(v1top(1:end-1),:),V2(v2top(1:end-1),:)},resolution,0,'linear'); Ft1 = fliplr(Ft1); 
[Ft2,Vt2] = regionTriMesh3D({V2(v2top(1:end-1),:),V3(v3top(1:end-1),:)},resolution,0,'linear');
[Ft3,Vt3] = regionTriMesh3D({V3(v3top(1:end-1),:),V4(v4top(1:end-1),:)},resolution,0,'linear'); Ft3 = fliplr(Ft3);

[Fb1,Vb1] = regionTriMesh3D({V1(v1bot(1:end-1),:),V2(v2bot(1:end-1),:)},resolution,0,'linear'); Fb1 = fliplr(Fb1);
[Fb2,Vb2] = regionTriMesh3D({V2(v2bot(1:end-1),:),V3(v3bot(1:end-1),:)},resolution,0,'linear'); Fb2 = fliplr(Fb2);
[Fb3,Vb3] = regionTriMesh3D({V3(v3bot(1:end-1),:),V4(v4bot(1:end-1),:)},resolution,0,'linear'); Fb3 = fliplr(Fb3);

if debug_plots
    figure
    hold on

    gpatch(F1,V1,'r');
    plot3V(V1(v1top,:),'r*')
    plot3V(V1(v1bot,:),'b*')

    gpatch(F2,V2,'b');
    plot3V(V2(v2top,:),'r*')
    plot3V(V2(v2bot,:),'b*')

    gpatch(F3,V3,'g');
    plot3V(V3(v3top,:),'r*')
    plot3V(V3(v3bot,:),'b*')

    gpatch(F4,V4,'m');
    plot3V(V4(v4top,:),'r*')
    plot3V(V4(v4bot,:),'b*')

    gpatch(Ft1,Vt1,'g');
    gpatch(Ft2,Vt2,'g');
    gpatch(Ft3,Vt3,'g');
    gpatch(Fb1,Vb1,'g');
    gpatch(Fb2,Vb2,'g');
    gpatch(Fb3,Vb3,'g');
    
    % axisGeom
end


%%
Ff = {F1,F2,F3,F4,Ft1,Ft2,Ft3,Fb1,Fb2,Fb3};
Vv = {V1,V2,V3,V4,Vt1,Vt2,Vt3,Vb1,Vb2,Vb3};
[F,V,~] = joinElementSets(Ff,Vv);
[F,V] = patchCleanUnused(F,V);
[F,V] = mergeVertices(F,V);


%% Create ellipsoid plaque

if percent_stenosis > 0
    rx_outer = bump_length * 0.22;
    ry_outer = percent_stenosis*lumen_radius*0.9;
    rz_outer = lumen_radius * 0.55;
    
    rx_inner = rx_outer - plq_inner_thickness;
    ry_inner = ry_outer - plq_inner_thickness;
    rz_inner = rz_outer - plq_inner_thickness;
    
    cx = buffer_length + (region_length/2);
    cy = -(lumen_radius - percent_stenosis*lumen_radius)/2 - 1.5*media_thick;
    cz = 0;
    
    
    % n is mesh resolution
    n = 30; 
    [X_outer, Y_outer, Z_outer] = ellipsoid(cx, cy, cz, rx_outer, ry_outer, rz_outer, n);
    [X_inner, Y_inner, Z_inner] = ellipsoid(cx, cy, cz, rx_inner, ry_inner, rz_inner, n);
    
    Outer_p = surf2patch(X_outer, Y_outer, Z_outer, 'triangles');
    Inner_p = surf2patch(X_inner, Y_inner, Z_inner, 'triangles');
    Fp_outer = Outer_p.faces; Vp_outer = Outer_p.vertices;
    Fp_inner = Inner_p.faces; Vp_inner = Inner_p.vertices;
    
    Vp_outer = applyCurvature(Vp_outer,kappa,curve_dir);
    Vp_inner = applyCurvature(Vp_inner,kappa,curve_dir);
    
    opt = struct();
    opt.pointSpacing=resolution;
    [Fp_outer, Vp_outer] = ggremesh(Fp_outer,Vp_outer,opt);
    [Fp_inner, Vp_inner] = ggremesh(Fp_inner,Vp_inner,opt);
    
    
    % Combine plaque mesh with everything else
    [Ft,Vt,Ct] = joinElementSets({F,Fp_outer,Fp_inner},{V,Vp_outer,Vp_inner});
else
    [Ft,Vt,Ct] = joinElementSets({F},{V});
end

%%
% Tetrahedralize
v_shift = [min(Vt(:,1)) 0 0];

if percent_stenosis > 0
    inner_region_V = mean(Vp_inner,1)-v_shift;
    [~,ind_outer] = min(sqrt(sum((Vp_inner(1,:)-Vp_outer).^2,2)));
    outer_region_V = mean([Vp_inner(1,:); Vp_outer(ind_outer,:)])-v_shift;
    [~,ind_other] = min(sqrt(sum((Vp_outer(ind_outer,:)-V4).^2,2)));
    other_region_V = mean([Vp_outer(ind_outer,:);V4(ind_other,:)])-v_shift;
end

inputStruct.Faces = Ft;
inputStruct.Nodes = Vt-v_shift;
inputStruct.faceBoundaryMarker = Ct;

if percent_stenosis > 0
    inputStruct.regionPoints = [inner_region_V;outer_region_V;other_region_V];
    inputStruct.regionA=[0.5,0.5,0.5];
else
    inputStruct.regionPoints = [];
    inputStruct.regionA = [];
end

inputStruct.holePoints = [];
inputStruct.stringOpt = sprintf('-pq%f/%fYA',1.05,0.0);
inputStruct.tetType   = 'tet10'; %Set desired element type
[mesh] = runTetGen(inputStruct);
inds = tesBoundary(mesh.faces);
mesh.facesBoundary = mesh.faces(inds,:);

u1 = p1(2,:) - mean(p1,1);   u2 = p2(2,:) - mean(p2,1);
v1 = p1(3,:) - mean(p1,1);   v2 = p2(3,:) - mean(p2,1);
w1 = cross(u1,v1);           w2 = cross(u2,v2);
w1 = w1/norm(w1);            w2 = w2/norm(w2);

mesh.surfs = separateMeshSurfacesNew(mesh.facesBoundary,mesh.nodes,[mean(V1(v1top,:),1);mean(V1(v1bot,:),1)],[-w1;-w2],resolution); % issue resolved for 7.04 sims, delete the resolution/2
mesh.ec = getCentroids(mesh.nodes,mesh.elements);
mesh.element_volume = getElementVolumes(mesh.elements,mesh.nodes);

% Stenosis flag for writing to FEB file
if percent_stenosis > 0
    mesh.stenosis = 1;
else
    mesh.stenosis = 0;
end

%%
% Assign materials
eid = mesh.elementMaterialID;
mesh.elementMaterialID = zeros(size(eid));

if percent_stenosis > 0
    [~,inner_plq_ind] = min(sqrt(sum((inner_region_V-mesh.ec).^2,2)));
    mesh.elementMaterialID(eid==eid(inner_plq_ind)) = 4; % inner plaque

    outer_plq_ind = find(mesh.elementMaterialID~=4 & sum(ismember(mesh.elements,unique(mesh.elements(mesh.elementMaterialID==4,:))),2)>=6,1);
    mesh.elementMaterialID(eid==eid(outer_plq_ind)) = 5; % outer plaque
end

[~,ind] = min(sqrt(sum((V1(v1top(1),:)-V2).^2,2)));
[~,iind] = min(sqrt(sum((mean([V2(ind,:);V1(v1top(1),:)],1)-mesh.ec).^2,2)));
mesh.elementMaterialID(eid==eid(iind)) = 1; % Media

[~,ind] = min(sqrt(sum((V2(v2top(1),:)-V3).^2,2)));
[~,iind] = min(sqrt(sum((mean([V3(ind,:);V2(v2top(1),:)],1)-mesh.ec).^2,2)));
mesh.elementMaterialID(eid==eid(iind)) = 2; % Adventitia

[~,ind] = min(sqrt(sum((V3(v3top(1),:)-V4).^2,2)));
[~,iind] = min(sqrt(sum((mean([V4(ind,:);V3(v3top(1),:)],1)-mesh.ec).^2,2)));
mesh.elementMaterialID(eid==eid(iind)) = 3; % Sleeve

end
