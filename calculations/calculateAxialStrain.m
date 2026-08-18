% Calculate axial tensile and compressive strain
% Inputs: strain_path --> .csv file of strain values output from FEBio sim
%         mesh --> mesh generation parameters & data
%         opts --> simulation parameters
% Output: axial strain (tensile and compressive) values per element for each time step

function [axialStrain] = calculateAxialStrain(strain_path, mesh, opts)

if isnumeric(strain_path)
    data = strain_path;
else
    data = readmatrix(strain_path);
end
data(isnan(data(:,1)),:) = [];

elem_id = data(:,1);
exx = data(:,2);
eyy = data(:,3);
ezz = data(:,4);
exy = data(:,5);
eyz = data(:,6);
exz = data(:,7);

% Calculate element centroids
elem_centroids = getCentroids(mesh.nodes, mesh.elements);

% Build centerline
npts = 300;
x = linspace(min(mesh.nodes(:,1)), max(mesh.nodes(:,1)), npts)';
C0 = [x, zeros(npts,2)];
C  = applyCurvature(C0, opts.kappa, opts.curve_dir);

% Create tangent vector field (gradient)
t = gradient(C);
t = t ./ vecnorm(t,2,2);

% Force consistent axial direction (+x)
flipMask = t(:,1) < 0;
t(flipMask,:) = -t(flipMask,:);

% Find closest centroid to centerline
Nc = size(C,1);
Ne = size(elem_centroids,1);
distMat = zeros(Ne,Nc);
for k = 1:Nc
    distMat(:,k) = vecnorm(elem_centroids - C(k,:),2,2);
end
[~, idx] = min(distMat,[],2);

% Pre-allocate
N = numel(elem_id);
t_ax   = zeros(N,1);
c_ax   = zeros(N,1);

% Loop through each element-timestep combination
for i = 1:N
    ei = elem_id(i);
    ti = t(idx(ei),:);

    % Assemble strain tensor
    E = [ exx(i) exy(i) exz(i);
          exy(i) eyy(i) eyz(i);
          exz(i) eyz(i) ezz(i) ];

    % Axial strain
    ea = ti * E * ti'; % Apply transformation
    t_ax(i) = max(ea,0);
    c_ax(i) = min(ea,0);

end

% Reshape by timestep
M = unique(elem_id);
num_steps = numel(elem_id) / numel(M);

t_ax   = reshape(t_ax,   num_steps, [])';
c_ax   = reshape(c_ax,   num_steps, [])';
axialStrain  = struct('tensileStrain', t_ax,'compressiveStrain', c_ax);

end