function R = rotation_matrix(theta, phi)

v = 1-cos(theta);
R = [1-v* cos(phi)^2, -cos(phi)*sin(phi)* v, cos(phi)* sin(theta);
    -cos(phi)*sin(phi)*v, 1- v*sin(phi)^2, sin(phi)*sin(theta);
    -cos(phi)*sin(theta), -sin(phi)* sin(theta), cos(theta)];
end
