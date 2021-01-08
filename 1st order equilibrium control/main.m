clear
clc
rng(1)
global legend_flag;
global noise_sigma;


legend_flag = [];

ENV = env(true);        % generate envrionment
polyCell = cell_decom(ENV.V, 0.5);                     % cell decomposition
DYN = dyn();        % first order dynamic
noise_sigma = eye(size(DYN.A)) * 0.1;
for i = 1:length(polyCell)  % loop for all cell
    polyK{i} = find_controller(DYN.A,DYN.B,polyCell{i});  % find controller for every cell
end

p0 = [-1,2.5];  % init position
x0 = [p0(1), p0(2)]';

[t, x, u] = simulation(DYN.A,DYN.B,x0,polyK,polyCell, true, 'b', 200);             % simulation
