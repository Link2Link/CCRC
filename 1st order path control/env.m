function polyCell = env()
%% This file generates the 2D polygon envrionment for path control
global legend_flag;
n_p = 10; 
% rng(123);
rand_point = rand(n_p, 2)*10 - 5;
x = rand_point(:,1);
y = rand_point(:,2);
dt=DelaunayTri(x,y);
k=convexHull(dt);
p_v = [x(k),y(k)];
p = p_v(1:end-1,:);
p_v = p_v - mean(p);
p = p - mean(p);

polyCell = cell_decom(p, 0.5);
polyCell = polyCell(1:6);

