function plotPolygon(vertex)
vertex = [vertex;vertex(1,:)];
hold on
p1 = plot(vertex(:,1),vertex(:,2), 'b-', 'DisplayName','solving zone');
hold off
end

