% TimZ, MSc, 2021, Imperial College.
% 2/1/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw a circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs
% r (2x1 Integer) = Center of the circle
% d (Integer) = Radius of the circle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs
% None
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function draw_circle(r, d)
phi = linspace(0,2*pi);
x = d * cos(phi) + r(1);
y = d * sin(phi) + r(2);
plot(x,y,'Color',[0.5 0.5 0.5])
end