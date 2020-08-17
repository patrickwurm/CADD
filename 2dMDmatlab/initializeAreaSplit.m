function area_split = initializeAreaSplit(pos, npart, a, xbox, ybox, example_name, dx0, dy0)
%This function splits the tensile test piece into npart equal parts in
%length direction
%This splitting is later used to calculate the current area of the tensile
%test piece by summing over the area of the npart quadliterals
%The cutout in the middle of the test piece is accustomed for (see plot at the end of the function)
%Note that the function does not work for too large values of npart because
%of the cutout in the middle
%The condition is: xbox/npart >= cutout/2

global tol

if strcmp(example_name,'Tensile_Test')

ybox = 29*dy0;

dx = (xbox-12*dx0)/npart;
xl=6*dx0;
xr=xl+dx;

yb = 10*dy0;
yt = 19*dy0;

for i = 1:npart
        for j = 1:4
            if j==1
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xl)<tol & abs(pos(:,2)-yb)<tol));
            elseif j==2
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xl)<tol & abs(pos(:,2)-yt)<tol));
            elseif j==3
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xr)<tol & abs(pos(:,2)-yt)<tol));
            else
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xr)<tol & abs(pos(:,2)-yb)<tol));
            end
        end
xl = xl+dx;
xr = xr+dx;
end

elseif strcmp(example_name, 'Tensile_Test_rotated')

ybox = 17*dy0;

dx = (xbox-30*dx0)/npart;
xl=15*dx0;
xr=xl+dx;

yb = 6*dy0;
yt = 11*dy0;

for i = 1:npart
        for j = 1:4
            if j==1
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xl)<tol & abs(pos(:,2)-yb)<tol));
            elseif j==2
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xl)<tol & abs(pos(:,2)-yt)<tol));
            elseif j==3
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xr)<tol & abs(pos(:,2)-yt)<tol));
            else
                area_split.element{i}.atomlist(j) = find((abs(pos(:,1)-xr)<tol & abs(pos(:,2)-yb)<tol));
            end
        end
xl = xl+dx;
xr = xr+dx;

end
    
    
end

%For testing purposes:
if 0>1
    figure(55)
    scatter(pos(area_split.element{1}.atomlist,1),pos(area_split.element{1}.atomlist,2))
    axis equal
    hold on
    for i = 2:length(area_split.element)
        scatter(pos(area_split.element{i}.atomlist,1),pos(area_split.element{i}.atomlist,2))
    end
end

end



