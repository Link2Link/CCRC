function draw_set(DYN,ENV,K)
%% draw the invariant set
global legend_flag;
d = 0.1;
Ah = ENV.M{1};
bh = ENV.M{2};
global legend_flag;

[x,y] = meshgrid(-5:d:5,-5:d:5);
X = [x(:), y(:)]';
idx = min(Ah*X+bh)>=0;
X = X(:, idx);

idx = zeros(1, size(X,2));

h = waitbar(0,'calculating invariant set, wait...');
for i = 1:size(X,2)
    x0 = X(:, i);
    x0  = [x0(1), x0(2)]';
    [t, x] = simulation(DYN,K, x0, false, 'b', 10);             % simulation
    if min(min(Ah*x+bh)) >= 0
        idx(i) = 1;
    end
    waitbar(i / size(X,2))
end
close(h)

X = X(:, idx==1);
x = X(1,:)';
y = X(2,:)';
dt=DelaunayTri(x,y);
k=convexHull(dt);
p_v = [x(k),y(k)];
hold on
h1 = plot(p_v(:, 1), p_v(:, 2), 'g', 'DisplayName','invariant set');
legend_flag = [legend_flag, h1];
hold off
end

