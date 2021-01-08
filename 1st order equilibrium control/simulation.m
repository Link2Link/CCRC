function [t, x, u] = simulation(A,B,x0,polyK,polyCell, drawflag, color, T)
global legend_flag;
global noise_sigma;
sigma = noise_sigma;

num = length(polyCell);
% T = 100;
dt = 0.1;
x = x0;
x_log = x0;
t_log = 0;
u_log = zeros(size(polyK{1}, 1), 1);
for t = 0:dt:T
    % decided poly
    for i = 1:num
       Ah = polyCell{i}.Ah;
       bh = polyCell{i}.bh;
       flag = sign(min(Ah*x+bh));
       if flag >= 0
           K = polyK{i};
           break
       end
    end
    u = K*(x + sigma * randn(2,1));
    [ts, x] = ode45(@(ts,x) A*x+B*u, [t,t+dt], x);
    x = x(end, :)';
    x_log = [x_log, x];
    t_log = [t_log, ts(end)];
    u_log = [u_log, u];
end

x = x_log;
t = t_log;
u = u_log;
if drawflag
    hold on
    h1 = plot(x(1,:),x(2,:), color, 'DisplayName','running trace');
    h2 = plot(x(1,1),x(2,1), ['*',color], 'DisplayName','start point');
    hold off
    legend_flag = [legend_flag, h1,h2];
end