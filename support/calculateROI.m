% Find elements in region of interest
% Inputs: mesh --> mesh generation parameters & data
%         opts --> simulation parameters
% Outputs: maskROI --> logical mask of elements in/out of ROI

function maskROI = calculateROI(mesh, opts)

p1_curve = applyCurvature([opts.buffer_length-1,0,0;opts.buffer_length,0,0;opts.buffer_length+1,0,0],opts.kappa,opts.curve_dir);
v = p1_curve(3,:) - p1_curve(1,:);
v = v./norm(v);

p2_curve = applyCurvature([(opts.buffer_length+opts.region_length)-1,0,0;opts.buffer_length+opts.region_length,0,0;(opts.buffer_length+opts.region_length)+1,0,0],opts.kappa,opts.curve_dir);
u = p2_curve(3,:) - p2_curve(1,:);
u = u./norm(u);

centroids = getCentroids(mesh.nodes,mesh.elements);

[~, side1] = projectPointOntoPlane(centroids, v, p1_curve(2,:));
[~, side2] = projectPointOntoPlane(centroids, u, p2_curve(2,:));
    
maskROI = side1>=0 & side2<=0;

end