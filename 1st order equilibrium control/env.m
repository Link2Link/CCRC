function env_struct = env(draw_flag)
%ENV create the polygon envrionment
%   input: num of vertex
%   output: envrionment struct

global legend_flag;

%% generate N random points
N = 10;
rand_point = rand(N, 2)*10 - 5;
x = rand_point(:,1);
y = rand_point(:,2);

%% calc the convex env from the random points
dt = delaunayTriangulation(x,y);
k=convexHull(dt);
p_v = [x(k),y(k)];          % p_v contains vertex+1 points
p = p_v(1:end-1,:);         % p containts vertex

p_v = p_v - mean(p);
p = p - mean(p);

%% tranfer vertex to matrix
[Ah, bh] = vex2h(p);
% Ah = [Ah, zeros(size(Ah,1), 2)];

%% data struct
env_struct = struct('Name', 'polygon');
env_struct.V = p;
env_struct.V_ = p_v;
env_struct.M = {Ah, bh};

if draw_flag
    %% plot env
    x_max = max(p(:,1));
    x_min = min(p(:,1));
    y_max = max(p(:,2));
    y_min = min(p(:,2));

    p1 = plot(p_v(:,1),p_v(:,2), 'o', 'DisplayName','vertices');
    hold on
    p2 = plot(p_v(:,1),p_v(:,2), 'r', 'DisplayName','polygon bound');
    axis([x_min-1 x_max+1 y_min-1 y_max+1])
    legend_flag = [legend_flag, p1, p2];

    xlabel('position $x$','Interpreter','LaTex');
    ylabel('position $y$','Interpreter','LaTex')
    hold off
end

end

