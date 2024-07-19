% TimZ, MSc, 2021, Imperial College.
% 9/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw a line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% r (2x1 Integer) = The centre of the line
% DOA (Integer) = Angle between the line and the positive x axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_line(r, DOA)

% Get the points on the lines
x = linspace(-1000,1000);
y = zeros(size(x));
points = [x; y];

% Rotation
rotation = [cos(DOA) -sin(DOA); sin(DOA) cos(DOA)];  % Rotation matrix
points = rotation * points;

% Shift
points = points + r;  % Shift to the centre

plot(points(1,:), points(2,:), 'Color',[0.5 0.5 0.5]);
end