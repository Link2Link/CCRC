function d = dis(A, b, x)
% calc the distance from given x to given line Ax+b = 0
% calc distance
d = A*x + b;
d = d / norm(A);

end

