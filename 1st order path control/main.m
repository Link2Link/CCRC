clear
clc
rng(1)

global legend_flag;
legend_flag = [];

global noise_sigma;
noise_sigma = 1;

fig = figure();

%% 
ENV = env(true);        % generate envrionment
DYN = dyn();        % first order dynamic
for i = 1:length(ENV.sub_polygon)                   % loop for all cell
    [Control{i}.K, Control{i}.D] = find_controller(DYN, ENV.sub_polygon(i));  % find controller for every cell
end

%% simulation 
p0 = [-1,2.5];  % init position
x0 = [p0(1), p0(2)]';

global git_flag
git_flag = 1;
x0 = ENV.sub_polygon(1).sample();
simulation(DYN, ENV,[1;3], Control, true, 'b', 50);             % simulation
