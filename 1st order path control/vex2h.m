function [Ah, bh] = vex2h(p)
%% fransfer from vetex p to poly function Ah bh
n_h = size(p, 1);
p_mid = mean(p);
p = [p;p(1,:)];
% cvx_begin quiet
%     variable Ah(n_h,2)
%     minimize(0)
%     subject to
%         for i = 1:n_h
%             x1 = p(i, :);
%             x2 = p(i+1, :);
%             Ah(i,:)*x1' + 1 == 0;
%             Ah(i,:)*x2' + 1 == 0;
%         end
% cvx_end
Ah = zeros(n_h, 2);
bh = zeros(n_h, 1);
for i = 1:n_h
    [A,b] = point2line(p(i, :), p(i+1, :));
    Ah(i, :) = A;
    bh(i, :) = b;
end


% for i = 1:size(Ah,1)
%     bh(i) = bh(i) / norm(Ah(i,:));
%     Ah(i,:) = Ah(i,:) / norm(Ah(i,:));
% end
s = sign(Ah*p_mid' + bh);
Ah = Ah.*s;
bh = bh.*s;
end