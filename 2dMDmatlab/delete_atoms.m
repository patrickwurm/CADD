function [pos, padpos] = delete_atoms(a, axmin, aymin, axmax, aymax, pos, padpos, example_name, dx, dy)

global tol

if strcmp(example_name,'Tensile_Test')
    
    r = 12*dx/2;
    r1 = 0.95*r;
    
    xc1 = r;
    xc2 = xc1;
    xc3 = axmax-r;
    xc4 = xc3;
    
    yc1 = 0;
    yc4 = yc1;
    yc2 = aymax;
    yc3 = yc2;
    
    n_delete = 0;
    
    i=1;
    while i <= length(pos(:,1))
        x=pos(i,1);
        y=pos(i,2);
        deleted = 0;
        if x >= r && x <= 112*dx/2
            if y > 19*dy+tol || y < 10*dy-tol
                pos(i,:) = [];
                deleted = 1;
                n_delete = n_delete + 1;
            end
        elseif (x-xc1)^2+(y-yc1)^2 < r1^2 || (x-xc2)^2+(y-yc2)^2 < r1^2 || (x-xc3)^2+(y-yc3)^2 < r1^2 || (x-xc4)^2+(y-yc4)^2 < r1^2
            pos(i,:) = [];
            deleted = 1;
            n_delete = n_delete + 1;
        end
        if deleted == 1
            %do nothing
        else
            i=i+1;
        end
    end
    
    disp(['-> I deleted ',num2str(n_delete),' atoms. There are now: ',num2str(length(pos(:,1))),' atoms.'])
    
elseif strcmp(example_name,'Tensile_Test_rotated')

    r = 10*dx;
    r1 = 0.95*r;
    
    xc1 = r;
    xc2 = xc1;
    xc3 = axmax-r;
    xc4 = xc3;
    
    yc1 = 0;
    yc4 = yc1;
    yc2 = aymax;
    yc3 = yc2;
    
    n_delete = 0;
    
    i=1;
    while i <= length(pos(:,1))
        x=pos(i,1);
        y=pos(i,2);
        deleted = 0;
        if x >= r && x <= axmax-r
            if y > 11*dy+tol || y < 6*dy-tol
                pos(i,:) = [];
                deleted = 1;
                n_delete = n_delete + 1;
            end
        elseif (x-xc1)^2+(y-yc1)^2 < r1^2 || (x-xc2)^2+(y-yc2)^2 < r1^2 || (x-xc3)^2+(y-yc3)^2 < r1^2 || (x-xc4)^2+(y-yc4)^2 < r1^2
            pos(i,:) = [];
            deleted = 1;
            n_delete = n_delete + 1;
        end
        if deleted == 1
            %do nothing
        else
            i=i+1;
        end
    end
    
    disp(['-> I deleted ',num2str(n_delete),' atoms. There are now: ',num2str(length(pos(:,1))),' atoms.'])
    
end

end

