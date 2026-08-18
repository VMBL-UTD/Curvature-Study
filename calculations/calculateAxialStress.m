% Calculate axial tensile and compressive stress
% Inputs: stress_path --> .csv file of stress values output from FEBio sim
%         mesh --> mesh generation parameters & data
%         opts --> simulation parameters
% Output: axial stress (tensile and compressive) values per element for each time step

function [axialData] = calculateAxialStress(stress_path, mesh, opts)

if isnumeric(stress_path)
    data = stress_path;
else
    data = readmatrix(stress_path);
end
data(isnan(data(:,1)),:) = [];

elem_id = data(:,1);
sx  = data(:,2);
sy  = data(:,3);
sz  = data(:,4);
sxy = data(:,5);
syz = data(:,6);
sxz = data(:,7);

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

    % Assemble stress tensor
    S = [ sx(i)  sxy(i) sxz(i);
          sxy(i) sy(i)  syz(i);
          sxz(i) syz(i) sz(i) ];

    % Axial stress
    sa = ti * S * ti'; % Apply transformation
    t_ax(i) = max(sa,0);
    c_ax(i) = min(sa,0);

end

% Reshape by timestep
M = unique(elem_id);
num_steps = numel(elem_id) / numel(M);

t_ax   = reshape(t_ax,   num_steps, [])';
c_ax   = reshape(c_ax,   num_steps, [])';
axialData  = struct('tensileStress', t_ax, 'compressiveStress', c_ax);

end