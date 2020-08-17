function [ in_triangle ] = triangleTest( nodes, position_PA )
%Test if a point with coordinates "position_PA" lies in the triangle 
%formed by the 3 nodes "nodes".

a = [nodes(1).material_coordinates, 0];
b = [nodes(2).material_coordinates, 0];
c = [nodes(3).material_coordinates, 0];

global tol
criteria = tol*1e-35; %the tolerance here must be very tight to ensure proper association of the pad atoms

criteria = tol*1e-35;

p = [position_PA, 0];

% activate this for a visualization of the pad atoms and the elements to
% be tested
if 0>1
    figure(71)
    line([a(1),b(1)], [a(2), b(2)], 'Color', 'r');
    hold on;
    line([b(1),c(1)], [b(2), c(2)], 'Color', 'r');
    line([a(1),c(1)], [a(2), c(2)], 'Color', 'r');
    plot(p(1),p(2),'k^');
    hold off
    drawnow
end

if sameSide(p, a, b, c) && sameSide(p, b, a, c) && sameSide (p, c, a, b)
    in_triangle = true;
else
    in_triangle = false;
end

    function on_same_side = sameSide(p1, p2, a, b)
        cp1 = cross(b-a, p1-a);
        cp2 = cross(b-a, p2-a);
        if dot(cp1,cp2) >= -criteria
            on_same_side = true;
        else
            on_same_side = false;
        end
    end

end

