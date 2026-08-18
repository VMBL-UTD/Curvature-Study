% Bend straight coordinates along specifed curvature direction
% Inputs: V --> Nx3 matrix of vertex coordinates [x, y, z]
%         R --> curvature magnitude, kappa
%         curvatureVec --> 1x3 vector specifying direction of curvature
% Outputs: V_curved --> Nx3 matrix of vertex coordinates after transform

function [V_curved] = applyCurvature(V, R, curvatureVec)

if R==0
    V_curved = V;
    return
else
    R = 1 / (R * 0.001);
end

% Convert curvature direction into unit vector
curvDir = curvatureVec(:)' / norm(curvatureVec);

% Set tube's axis of extrusion
tubeAxis = [1 0 0];

w = cross(tubeAxis, curvDir);
if dot(cross(curvDir, w), tubeAxis) < 0
    w = -w;
end
w = w / norm(w);

% Define coordinate space within plane of curvature (orthonormal basis)
x_local = V * curvDir';                     % treat direction of curvature as x
y_local = V * w';                           % y is axis orthogonal to both tube 
z_local = V * tubeAxis';                    % z is the direction along artery 

% Define new coordinates relative to plane of curvature
theta = z_local / R;                        % calculate angle of arc (arc_length/R)
x_new = (R + x_local) .* cos(theta) - R;    % Subtracting R to keep tube at origin
z_new = (R + x_local) .* sin(theta);

V_curved = x_new .* curvDir + y_local .* w + z_new .* tubeAxis;

end