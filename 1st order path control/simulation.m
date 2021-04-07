function [t, x, u] = simulation(DYN,ENV, x0,control, drawflag, color, T)
global legend_flag;
global noise_sigma;
global git_flag
sigma = noise_sigma;

num = length(ENV.sub_polygon);
A = DYN.A;
B = DYN.B;
dt = 0.1;
x = x0;
x_log = x0;
t_log = 0;
u_log = zeros(2, 1);

hold on
state = plot(x0(1), x0(2), 'ro');
trajctory = animatedline;
    addpoints(trajctory,x(1),x(2));
h2 = plot(x(1,1),x(2,1), ['*',color], 'DisplayName','start point');
for t = 0:dt:T
    % decided poly
    for i = 1:num
        if ENV.sub_polygon(i).in_polygon(x)
            K = control{i}.K;
            D = control{i}.D;
            break
        end
    end
    
    % apply control
    u = K*(x + sigma * randn(2,1))+ D;
    [ts, x] = ode45(@(ts,x) A*x+B*u, [t,t+dt], x);
    x = x(end, :)';
    state.XData = x(1);
    state.YData = x(2);
    addpoints(trajctory,x(1),x(2));
    x_log = [x_log, x];
    t_log = [t_log, ts(end)];
    u_log = [u_log, u];
    drawnow
    frame = getframe;
    [imind,cm] = rgb2ind(frame.cdata,256);
    if git_flag
        imwrite(imind,cm,'rename.gif','gif','LoopCount',Inf,'DelayTime',dt/10000);
        git_flag = 0;
    else
        imwrite(imind,cm,'rename.gif','gif','WriteMode','append' ,'DelayTime',dt/10000);
    end
%     pause(dt/10);
end

x = x_log;
t = t_log;
u = u_log;
