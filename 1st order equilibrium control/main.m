clear
clc
rng(1)

global legend_flag;
legend_flag = [];

global noise_sigma;
noise_sigma = 0.5;

fig = figure();

%% 
ENV = env(true);        % generate envrionment
DYN = dyn();        % first order dynamic
for i = 1:length(ENV.sub_polygon)                   % loop for all cell
    Control{i} = find_controller2(DYN, ENV.sub_polygon(i));  % find controller for every cell
end

%% simulation 
p0 = [-1,2.5];  % init position
x0 = [p0(1), p0(2)]';

global git_flag
git_flag = 1;
for i = 1:length(ENV.sub_polygon)
    x0 = ENV.sub_polygon(i).sample();
    simulation(DYN, ENV,x0, Control, true, 'b', 10);             % simulation
end
