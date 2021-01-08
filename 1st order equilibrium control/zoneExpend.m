function SolvingZone = zoneExpend(SolvingZone,ENV)
N = 1e5;
xq = rand(N, 1)*10 - 5;
yq = rand(N, 1)*10 - 5;
xv = ENV.V_(:, 1);
yv = ENV.V_(:, 2);

x = SolvingZone(:, 1);
y = SolvingZone(:, 2);
dt = delaunayTriangulation(x,y);
k=convexHull(dt);
k = k(end:-1:1);
p_v = [x(k),y(k)];          % p_v contains vertex+1 points
xv = [xv', NaN, x(k)'];
yv = [yv', NaN, y(k)'];

in = inpolygon(xq,yq,xv,yv);
if sum(in) > 0
    new_xs = xq(in);
    new_ys = yq(in);

    x = [x; new_xs(1)];
    y = [y; new_ys(1)];
    dt = delaunayTriangulation(x,y);
    k=convexHull(dt);
    p_v = [x(k),y(k)];          % p_v contains vertex+1 points
    p = p_v(1:end-1,:);         % p containts vertex

    SolvingZone = p;
end
end

