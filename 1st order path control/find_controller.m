function K = find_controller(A,B,polyCellStruct)
global noise_sigma;
r = 2;
Ah = polyCellStruct.Ahi_ref;
bh = polyCellStruct.bhi_ref;
z = polyCellStruct.Az_ref;
bz = polyCellStruct.bz_ref;
x_v = polyCellStruct.vertex_ref;


n_h = size(Ah, 1);
n_p = size(x_v, 1);
n_x = size(A, 2);

xk = [x_v(:, 1), x_v(:, 2)];

sigma = noise_sigma;
alpha_hi = ones(1,n_h);
beta = [1, 1];
t = 2;
ti = t;
tv = t;
gamma = 1;

eta = @(x) log(1/gamma * dis(Ah, bh, x) + 1); 
eta_v = 0.1;

if ~isempty(Ah)
    Ai = alpha_hi(1)*Ah + alpha_hi(2)*Ah*A;
end
Az = beta(1)*z + beta(2)*z;
alphai = alpha_hi(1) * bh;
beta_z = beta(1)* bz;

cvx_begin quiet
    variables K(2,2)
    variable relax(n_h)
    variable relax_v
    minimize(norm(K) + 1e3 * sum(relax) + 1e1 * relax_v)
    subject to
        %% chance CBP constraint
        for i = 1:n_h
            for k = 1:n_p
                x = xk(k, :)';
                Ki = Ah(i,:)*B*K;
                Ki_p = Ai(i,:) + Ki;
                Kik_p = Ki_p * x;
                eta_ = eta(x);
                Kik_p'*Kik_p + 2*(alphai(i) - ti)*Ki_p*x - ti^2*eta_(i)...
                    + Ki*sigma*Ki' + (alphai(i) - ti)^2 <= relax(i);
            end
        end
        %% chance CLF constraint
        for k = 1:n_p
            x = xk(k, :)';
            Kz = z*B*K;
            Kz_p = Az + Kz;
            Kzk_p = Kz_p * x;
            Kzk_p'*Kzk_p + 2*(beta_z - tv)*Kz_p*x - tv^2*eta_v...
                + Kz*sigma*Kz' + (beta_z + tv)^2 <= relax_v;
        end
        relax >= 0;
        relax_v >= 0;
cvx_end
relax
end