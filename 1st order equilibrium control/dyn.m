function  dyn_struct = dyn()
%% first order dynamic

A = zeros(2);
B = eye(2);

dyn_struct = struct('A',A, 'B',B);