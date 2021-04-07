function [K, D] = find_controller(DYN, ENV)
global noise_sigma;

r = 1;
A = DYN.A;
B = DYN.B;
Ah = ENV.h_ref{1};
bh = ENV.h_ref{2};
Az = ENV.z_ref{1};
bz = ENV.z_ref{2};
Ax = ENV.M_ref{1};
bx = ENV.M_ref{2};

xe = ENV.p_ref';
D = 0;


n_h = size(Ah, 1);
n_z = size(Az, 1);
n_x = size(A, 2);
n_ax = size(Ax, 1);

sigma = noise_sigma;
alpha_hi = ones(1,n_h);
beta = 1;

P = eye(n_x);
C = eye(n_x);

cvx_begin quiet
    variables K(2,2)
    variable Q(2,2) symmetric
    variable lambda(n_ax, n_h)
    variable lambda_z(n_ax, n_z)
    variable punish_h(n_h)
    variable punish_z(n_z)
    minimize(norm(K))
    subject to
        %% chance CBP constraint
        for i = 1:n_h
            Lg = Ah(i, :)*B;
            T = Ah(i, :)*A + alpha_hi(i)*Ah(i, :);
            K1 = T + Lg*K*C;
            K2 = Lg*K;
            alpha_ = alpha_hi(i)*bh(i);
            Gamma = sigma;
            K2*Gamma*K2' - Lg*D - alpha_ <= -lambda(:, i)'*bx;
            K1 - lambda(:, i)'*Ax == 0;
            lambda(:, i)' >= 0;
        end

        - Q ==  ((A+B*K)'*P + P*(A+B*K) + beta*P);  
        Q >= 0;
cvx_end

D = -K*xe;

end


function d = dis(A, b, x)
% calc the distance from given x to given line Ax+b = 0
% calc distance
d = A*x + b;
d = d / norm(A);

end