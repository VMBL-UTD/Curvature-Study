function distance = pointToPlaneDistance(P1, P2, P3, Q)
    % Inputs:
    % P1, P2, P3 - Three points defining the plane (1x3 vectors)
    % Q - The point off the plane (1x3 vector)
    % Output:
    % distance - The distance between point Q and the plane
    
    % Calculate two vectors in the plane
    v1 = P2 - P1;
    v2 = P3 - P1;
    
    % Calculate the normal vector to the plane
    normal = cross(v1, v2);
    
    % Normalize the normal vector
    normal = normal / norm(normal);
    
    % Plane equation: ax + by + cz + d = 0
    % d is calculated as: -dot(normal, P1)
    d = -dot(normal, P1);
    
    % Calculate the distance from Q to the plane
    distance = abs(dot(normal, Q) + d) ./ norm(normal);
end