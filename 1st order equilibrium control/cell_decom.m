function polyCell = cell_decom(p, ratio)
%% cell decomposition
%% ratio is the ratio for cell decomposition
global legend_flag;


mid = mean(p);
p = [p;p(1,:)];
p_m = (p - mid)*ratio + mid;

for num = 1:size(p,1)-1
    points = [p(num, :);p(num+1, :);p_m(num+1, :);p_m(num, :)];
    polyCell{num}.p = points;
    [Ah, bh] = vex2h(points);
    polyCell{num}.Ah = Ah;
    polyCell{num}.bh = bh;
    idx = min(abs(Ah*points(1:2,:)' + bh) < 1e-3, [], 2);
    Ahi = Ah(idx, :);
    bhi = bh(idx, :);
    polyCell{num}.Ahi = Ahi;
    polyCell{num}.bhi = bhi;
end

polyCell{size(p,1)}.p = p_m(1:end-1,:);
[Ah, bh] = vex2h(p_m(1:end-1,:));
polyCell{size(p,1)}.Ah = Ah;
polyCell{size(p,1)}.bh = bh;
polyCell{size(p,1)}.Ahi = [];
polyCell{size(p,1)}.bhi = [];

for i = 1:length(polyCell)
    polyCell{i}.Ah = [polyCell{i}.Ah(:, 1), polyCell{i}.Ah(:, 2)];
    if ~isempty(polyCell{i}.Ahi)
        polyCell{i}.Ahi = [polyCell{i}.Ahi(1), polyCell{i}.Ahi(2)];
    end
end

hold on
h = plot(p_m(:,1), p_m(:,2), 'r--',  'DisplayName','cell decomposition');
for num = 1:size(p,1)-1
    plot(polyCell{num}.p(:,1), polyCell{num}.p(:,2), 'r--');
    s=sprintf('cell %d',num);
    text(mean(polyCell{num}.p(:,1)), mean(polyCell{num}.p(:,2)), s)
end
hold off
legend_flag = [legend_flag, h];

end