clear
close all
clc
rng(1)
global legend_flag;
global noise_sigma;
legend_flag = [];
polyCell = env;                                 % generate envrionment
[A,B] = dyn;
noise_sigma = eye(size(A)) * 0.01;
for i = 1:length(polyCell)  % loop for all cell
    polyK{i} = find_controller(A,B,polyCell{i});  % find controller for every cell
end
% polyK{4} = 4 * polyK{4};
rng(2)

p0 = [1, 3];  % init position
v0 = [0, 0];     % init velocity
x0 = [p0(1), p0(2)]';
[t, x, u] = simulation(A,B,x0,polyK,polyCell, true, 'b');             % simulation

% hold on;
% h = plot(p0(1), p0(2), 'b*', 'DisplayName','start point');
% legend_flag = [legend_flag, h];
legend(legend_flag)

xlabel('position x')
ylabel('position y')