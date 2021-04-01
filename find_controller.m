function K = find_controller(DYN, ENV)
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
alpha_hi = ones(1,n_h);
t = 2;
ti = t;
tv = t;
gamma = 1;

eta = @(x) log(1/gamma * dis(Ah, bh, x) + 1); 
eta_v = 0.1;
alphai = alpha_hi * bh;
beta = 1;
P = eye(n_x);

if ~isempty(Ah)
    Ai = alpha_hi(1)*Ah;
end



cvx_begin quiet
    variables K(2,2)
    variable Q(2,2) symmetric
    variable relax(n_h)
    variable relax_v
    variable a
    minimize(1e4 * sum(relax) + 1e1 * relax_v)
    subject to
    %% chance CBP constraint
    for i = 1:n_h
        Ai = alpha_hi(1)*Ah;
        for k = 1:n_p
            x = xk(k, :)';
            Ki = Ah(i, :)*B*K;
            Ki_p = Ai(i, :) + Ki;
            Kik_p = Ki_p * x;
            pro = ti^2*eta(x);
            Kik_p'*Kik_p + 2*(alphai(i) - ti)*Ki_p*x - pro(i)...
                    + Ki*sigma*Ki' + (alphai(i) - ti)^2 <= relax(i);
        end
    end
    %% chance CLF constraint
    for k = 1:n_p
        x = xk(k, :)';
        - Q ==  ((A+B*K)'*P + P*(A+B*K) + beta*P);  
        Kv = P*B*K;
        Qk = Q*x;
        Kvk = Kv'*x;
        Qk'*x*x'*Qk + 4*Kvk'*sigma*Kvk + tv^2*(1 - eta_v) <= relax_v;
    end
    relax_v >= 0;
    relax >= 0;
    Q >= 0;
cvx_end
relax;
end


function d = dis(A, b, x)
% calc the distance from given x to given line Ax+b = 0
% calc distance
d = A*x + b;
d = d / norm(A);

end