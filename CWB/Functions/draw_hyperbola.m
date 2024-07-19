% TimZ, MSc, 2021, Imperial College.
% 2/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw a hyperbola (Ellipse)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% r_1 (2x1 Integer) = 1st focus of the ellipse
% r_2 (2x1 Integer) = 2nd focus of the ellipse
% d (Integer) = length of the major axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_hyperbola(r_1, r_2, d)

a = d/2;  % semi-major axes
c = norm(r_1-r_2)/2;  % Linear eccentricity
b = sqrt(c^2-a^2);  % semi-minor axes

% Get the points on the hyperbola (Ellipse)
y = linspace(-100,100);
if d >= 0
    x = -a * sqrt(1 + y.^2 / b^2);
else
    x = a * sqrt(1 + y.^2 / b^2);
end
points = [x; y];

% Rotation
% Determine the direction and rotate the hyperbola
r_1_2 = r_1 - r_2;
phi = atan(r_1_2(2) / r_1_2(1));  % Rotation angle
rot = [cos(phi) -sin(phi); sin(phi) cos(phi)];  % Rotation matrix
points = rot * points;

% Shift
r = (r_1 + r_2) / 2;  % Use the focuses to find the centre
points = points + r; % Shift to the centre

plot(points(1,:),points(2,:),'Color',[0.5 0.5 0.5]);
end