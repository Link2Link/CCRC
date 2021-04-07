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

%% data struct
env_struct = struct('Name', 'polygon');
env_struct.V = p;
env_struct.V_ = p_v;
env_struct.M = {Ah, bh};
env_struct.in_polygon = @(x) inside_polygon(x, Ah, bh);

%% plot

if draw_flag
    hold on
    %% plot env
    x_max = max(p(:,1));
    x_min = min(p(:,1));
    y_max = max(p(:,2));
    y_min = min(p(:,2));

    p1 = plot(p_v(:,1),p_v(:,2), 'o', 'DisplayName','vertices');
    p2 = plot(p_v(:,1),p_v(:,2), 'r', 'DisplayName','polygon bound');
    axis([x_min-1 x_max+1 y_min-1 y_max+1])
    legend_flag = [legend_flag, p1, p2];

    xlabel('position $x$','Interpreter','LaTex');
    ylabel('position $y$','Interpreter','LaTex')
    hold off
end

%% decomposition
ratio = 0.5;
mid = mean(env_struct.V);
p = env_struct.V_;
p_m = (p - mid)*ratio + mid;

for num = 1:size(p,1)-1
    sub_struct = struct('Name', 'sub_polygon');
    points = [p(num, :);p(num+1, :);p_m(num+1, :);p_m(num, :)];
    sub_struct.V = points;
    [Ah, bh] = vex2h(points);
    sub_struct.M = {Ah, bh};
    sub_struct.in_polygon = @(x) inside_polygon(x, Ah, bh);
    sub_struct.sample = @() sample(Ah, bh);
    
    p_ref = (points(2,:) + points(3,:)) / 2 - Ah(2, :)/norm(Ah(2, :))*1;
    sub_struct.p_ref = p_ref;
    
    points = points - p_ref;
    sub_struct.vertex_ref = points;
    [Ah_ref, bh_ref] = vex2h(points);
    sub_struct.M_ref = {Ah_ref, bh_ref};
    
    Ahi = Ah([1,3], :);
    bhi = bh([1,3], :);
    Az = Ah(2, :);
    bz = bh(2, :);
    
    Ahi_ref = Ah_ref([1,3], :);
    bhi_ref = bh_ref([1,3], :);
    Az_ref = Ah_ref(2, :);
    bz_ref = bh_ref(2, :);
    
    sub_struct.h = {Ahi, bhi};
    sub_struct.h_ref = {Ahi_ref, bhi_ref};
    
    sub_struct.z = {Az, bz};
    sub_struct.z_ref = {Az_ref, bz_ref};
    
    env_struct.sub_polygon(num) = sub_struct;
end

% num = size(p,1);
% sub_struct = struct('Name', 'sub_polygon');
% sub_struct.V = p_m(1:end-1,:);
% [Ah, bh] = vex2h(sub_struct.V);
% sub_struct.M = {Ah, bh};
% sub_struct.in_polygon = @(x) inside_polygon(x, Ah, bh);
% sub_struct.sample = @() sample(Ah, bh);
% Ahi = [];
% bhi = [];
% sub_struct.h = {Ahi, bhi};    
% env_struct.sub_polygon(num) = sub_struct;

%% plot
hold on
h = plot(p_m(:,1), p_m(:,2), 'r-',  'DisplayName','cell decomposition');
for num = 1:size(p,1)-1
    p = env_struct.sub_polygon(num).V;
    plot(p(:,1), p(:,2), 'r--');
    s=sprintf('cell %d',num);
%     text(mean(p(:,1)), mean(p(:,2)), s);
end
hold off
legend_flag = [legend_flag, h];



end

function out = inside_polygon(x, Ah, bh)
flag = Ah*x + bh > 0;
out = sum(flag) == length(flag);
end

function x=sample(Ah, bh)
x = [];
while isempty(x)
    rand_point = rand(2, 1)*10 - 5;
    if inside_polygon(rand_point, Ah, bh)
        x = rand_point;
    end
    
end
end


