function [ positions ] = genatoms(a, xmin, ymin, xmax, ymax, dothecut, dx, dy )

%Has capabilities to produce a Volterra cut via dothecut = 1 (was not needed in the end)
%Uses old version of genatoms (which produces a lattice that is rotated by 90°)

addpath('2dMDmatlab')

positions = genatoms_original(a, ymin, xmin, ymax, xmax, dy, dx);

positions_x = positions(:,2);
positions_y = positions(:,1);

positions(:,1) = positions_x;
positions(:,2) = positions_y;

%scatter(positions(:,1), positions(:,2));

xcutpos = (xmax-xmin)/2;
ycutpos = (ymax-ymin)*9/10;

burger = a;

% make a copy of atoms to be added later
saveatoms = [];
for i = 1:length(positions(:,1))
    if positions(i,1) < xcutpos
        if positions(i,2) >= max(positions(:,2))- 1.1*1/2*a;
            saveatoms = [saveatoms; positions(i,:)];
        end
    end
end
saveatoms(:,2) = saveatoms(:,2) + burger/2;

if dothecut
% do the Volterra cut
for i = 1:length(positions(:,1))
    if positions(i,1) < xcutpos
        if positions(i,2) > ycutpos
            positions(i,2) = positions(i,2) - burger/2;
        end
    end
    if positions(i,1) > xcutpos
        if positions(i,2) > ycutpos
            positions(i,2) = positions(i,2) + burger/2;
        end
    end
end

%add the "missing" atoms
positions = [positions; saveatoms];

% figure(1)
% scatter(positions(:,1), positions(:,2));
% axis equal
% drawnow
end

end