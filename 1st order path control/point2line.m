function [A,b] = point2line(p1,p2)
%POINT2LINE 
x1 = p1(1);
y1 = p1(2);
x2 = p2(1);
y2 = p2(2);

if abs(x1-x2) < 1e-5
    a = 1;
    b = 0;
    c = -x1
else
    k = (y2 - y1) / (x2 - x1);
    a = k;
    b = -1;
    c = -k*x1+y1;
end
A = [a,b];
b = c;

b = b / norm(A);
A = A / norm(A);