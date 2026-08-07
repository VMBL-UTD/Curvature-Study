function [F,V,vtop,vbot] = createStenosedLumen(opts)

if isfield(opts,'region_length'); region_length = opts.region_length; else; error("Field 'region_length' missing from input struct"); end
if isfield(opts,'buffer_length'); buffer_length = opts.buffer_length; else; error("Field 'buffer_length' missing from input struct"); end
if isfield(opts,'resolution'); resolution = opts.resolution; else; error("Field 'resolution' missing from input struct"); end
if isfield(opts,'lumen_radius'); lumen_radius = opts.lumen_radius; else; error("Field 'lumen_radius' missing from input struct"); end
if isfield(opts,'percent_stenosis'); percent_stenosis = opts.percent_stenosis; else; error("Field 'percent_stenosis' missing from input struct"); end
if isfield(opts,'kappa'); kappa = opts.kappa; else; error("Field 'kappa' missing from input struct"); end
if isfield(opts,'curve_dir'); curve_dir = opts.curve_dir; else; error("Field 'curve_dir' missing from input struct"); end

n = 20;
th = linspace(0,2*pi,n);

% End cap 1 profile
[y,z] = pol2cart(th,lumen_radius);
p1 = [zeros(n,1), y', z'];

% Region start profile
[y,z] = pol2cart(th,lumen_radius);
p2 = [ones(n,1)*buffer_length, y', z'];
m2 = [buffer_length,0,0];
n2 = cross(p2(1,:)-m2,p2(2,:)-m2); n2 = n2./norm(n2);

% Mid stenosis profile
r = lumen_radius-(lumen_radius*percent_stenosis);
[y,z] = pol2cart(linspace(0,2*pi,n),r);
ax = buffer_length+region_length/2;
p3 = [ones(n,1)*ax, y' - (r-lumen_radius), z'];
m3 = [ax,-(r-lumen_radius),0];
n3 = cross(p3(1,:)-m3,p3(2,:)-m3); n3 = n3./norm(n3);

% Region of interest end profile
[y,z] = pol2cart(th,lumen_radius);
ax = region_length + buffer_length;
p4 = [ones(n,1)*ax, y', z'];

% Region end profile
[y,z] = pol2cart(th,lumen_radius);
ax = region_length + 2*buffer_length;
p5 = [ones(n,1)*ax, y', z'];

% Create the mesh
na = 10;
cPar.closeLoopOpt = 1;
cPar.patchType = 'tri_slash';
cPar.numSteps = na;

% Initial buffer section
[Fi,Vi] = polyLoftLinear(p1,p2,cPar);

% Region of interest 1
[Fr1,Vr1] = sweepProfileCustom(p2,p3,n2,n3,na);

% Region of interest 2 (just flip the first region of interest)
Fr2 = Fr1;
Vr2 = Vr1;
Vr2(:,1) = -Vr2(:,1) + 2*max(Vr1(:,1));

% End buffer section
[Ff,Vf] = polyLoftLinear(p4,p5,cPar);

% Combine and clean
[F,V] = joinElementSets({Fi,Fr1,Fr2,Ff},{Vi,Vr1,Vr2,Vf});
[F,V] = mergeVertices(F,V);
[F] = patchCleanUnused(F,V);

% Apply curvature
V = applyCurvature(V,kappa,curve_dir);

% Remesh to get uniform triangles
optstr.pre.max_hole_area = 10000;
optstr.pre.max_hole_edges = 0;
optstr.pointSpacing = resolution;
[F,V] = ggremesh(F,V,optstr);

% Get mesh edge boundaries
[E,~,~] = patchBoundary(F);
optionStruct.outputType='label';
g1=tesgroup(E,optionStruct);
ve1 = edgeListToCurve(E(g1==1,:));
ve2 = edgeListToCurve(E(g1==2,:));

d1 = sqrt(sum(mean(V(ve1,:),1).^2,2));
d2 = sqrt(sum(mean(V(ve2,:),1).^2,2));
if d1<d2
    vtop = ve1;
    vbot = ve2;
else
    vtop = ve2;
    vbot = ve1;
end

pd1 = zeros(length(vtop),1);
for i=1:length(vtop)
    pd1(i,:) = pointToPlaneDistance(mean(p1,1),p1(2,:),p1(3,:),V(vtop(i),:));
end

pd2 = zeros(length(vbot),1);
vend = applyCurvature(p5,kappa,curve_dir);
for i=1:length(vbot)
    pd2(i,:) = pointToPlaneDistance(mean(vend,1),vend(2,:),vend(3,:),V(vbot(i),:));
end

V(vtop,:) = V(vtop,:)-[1,0,0].*pd1;
V(vbot,:) = V(vbot,:)+cross(vend(1,:)-mean(vend,1),vend(2,:)-mean(vend,1)).*pd2;

end


function [F,V] = sweepProfileCustom(p1,p2,n1,n2,na)
vc = struct();
for i=1:size(p1,1)
    vc(i).v = sweepCurveBezier(p1(i,:),p2(i,:),n1,n2,na,0.5);
end
vc = vertcat(vc.v);

cPar.closeLoopOpt=1;
cPar.patchType='tri_slash';
cPar.numSteps=1;

tmpf = cell(na-1,1);
tmpv = cell(na-1,1);

for i=1:na-1
    v1 = vc((1:na:na*size(p1,1))+(i-1),:);
    v2 = vc((1:na:na*size(p1,1))+i,:);

    [tmpf{i},tmpv{i}] = polyLoftLinear(v1,v2,cPar);
end
[F,V] = joinElementSets(tmpf,tmpv);
end

