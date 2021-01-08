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
    
    p_ref = (points(2,:) + points(3,:)) / 2;
    polyCell{num}.p_ref = p_ref;
    points = points - p_ref;
    polyCell{num}.vertex_ref = points;
    [Ah_ref, bh_ref] = vex2h(points);
    polyCell{num}.Ah_ref = Ah_ref;
    polyCell{num}.bh_ref = bh_ref;
    
end

polyCell{size(p,1)}.p = p_m(1:end-1,:);
[Ah, bh] = vex2h(p_m(1:end-1,:));
polyCell{size(p,1)}.Ah = Ah;
polyCell{size(p,1)}.bh = bh;

for i = 1:length(polyCell)-1
    polyCell{i}.Ah = [polyCell{i}.Ah(:, 1), polyCell{i}.Ah(:, 2)];
    polyCell{i}.Ah_ref = [polyCell{i}.Ah_ref(:, 1), polyCell{i}.Ah_ref(:, 2)];
%     if ~isempty(polyCell{i}.Ahi)
%         polyCell{i}.Ahi = [polyCell{i}.Ahi(1), 0, polyCell{i}.Ahi(2), 0];
%     end
    polyCell{i}.Ahi = polyCell{i}.Ah([1,3], :);
    polyCell{i}.bhi = polyCell{i}.bh([1,3]);
    polyCell{i}.Az = polyCell{i}.Ah(2, :);
    polyCell{i}.bz = polyCell{i}.bh(2);
    
    polyCell{i}.Ahi_ref = polyCell{i}.Ah_ref([1,3], :);
    polyCell{i}.bhi_ref = polyCell{i}.bh_ref([1,3]);
    polyCell{i}.Az_ref = polyCell{i}.Ah_ref(2, :);
    polyCell{i}.bz_ref = polyCell{i}.bh_ref(2);
    
end

hold on
h = plot(p_m(:,1), p_m(:,2), 'r-',  'DisplayName','boundary');

for num = 1:size(p,1)-1
    plot(polyCell{num}.p([1,2],1), polyCell{num}.p([1,2],2), 'r-');
    s=sprintf('cell %d',num);
    text(mean(polyCell{num}.p(:,1)), mean(polyCell{num}.p(:,2)), s)
    plot(polyCell{num}.p([2,3],1), polyCell{num}.p([2,3],2), 'b--');
end
hold off
legend_flag = [legend_flag, h];

end