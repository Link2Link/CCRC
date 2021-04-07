function K = find_controller2(DYN, ENV)
global noise_sigma;

r = 1;
A = DYN.A;
B = DYN.B;
Ah = ENV.h{1};
bh = ENV.h{2};
x_v = ENV.V;
n_h = size(Ah, 1);
n_p = size(x_v, 1);
n_x = size(A, 2);

xk = x_v;
sigma = noise_sigma;
alpha = 1;
eta_v = 0.1;

if ~isempty(Ah)
    Lg = Ah*B;
    T = Ah*A + alpha*Ah;
    alpha_ = alpha*bh;
    d = max(Ah*x_v'+bh');
    gamma = noise_sigma/d;
    eta_0 = 1/d;
    Gamma = gamma/eta_0;
end
Ax = ENV.M{1};
bx = -ENV.M{2};

beta = 1;
P = eye(n_x);

C = eye(n_x);
D = 0;


cvx_begin quiet
    variables K(2,2)
    variable Q(2,2) symmetric
    variable lambda(4,1)
    minimize(norm(K))
    subject to
    %% chance CBP constraint
    if ~isempty(Ah)
        K1 = T+Lg*K*C;
        K2 = Lg*K;
        K2*Gamma*K2' - K2*D + lambda'*bx - alpha_ <= 0;
        lambda'*Ax + K1 == 0;
        lambda >= 0;
    end
    - Q ==  ((A+B*K)'*P + P*(A+B*K) + beta*P);  
    Q >= 0;
cvx_end
end


function d = dis(A, b, x)
% calc the distance from given x to given line Ax+b = 0
% calc distance
d = A*x + b;
d = d / norm(A);

end