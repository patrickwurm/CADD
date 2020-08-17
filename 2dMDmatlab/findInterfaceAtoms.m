function interface_atom_list = findInterfaceAtoms(pos, xmin, ymin, xmax, ymax, a, example)

global tol

iatoms = 0;

if strcmp(example,'Example1')
    
    for i=1:length(pos(:,1))
        if (abs(pos(i,1)-xmin)<tol && pos(i,2)>=ymin) || (abs(pos(i,1)-xmax)<tol && pos(i,2)>=ymin) ...
                || (abs(pos(i,2)-ymin)<tol && (pos(i,1) >= xmin && pos(i,1) <= xmax))
            iatoms = iatoms + 1;
            interface_atom_list(iatoms,1) = i;
        end
    end
    
elseif strcmp(example,'Example2') || strcmp(example,'Example2_rotated')
    for i=1:length(pos(:,1))
        if (abs(pos(i,1)-xmin)<tol && pos(i,2)>=ymin)
            iatoms = iatoms + 1;
            interface_atom_list(iatoms,1) = i;
        end
    end
    
elseif strcmp(example,'Dislocation')
    
    for i=1:length(pos(:,1))
        if (abs(pos(i,1)-xmin)<tol && pos(i,2)>=ymin) || (abs(pos(i,1)-xmax)<tol && pos(i,2)>=ymin) ...
                || (abs(pos(i,2)-ymin)<tol && (pos(i,1) >= xmin && pos(i,1) <= xmax))
            iatoms = iatoms + 1;
            interface_atom_list(iatoms,1) = i;
        end
    end
    
elseif strcmp(example,'Radial_Pulse')
    
    for i=1:length(pos(:,1))
        if (abs(pos(i,1)-xmin)<tol && pos(i,2)>=ymin) || (abs(pos(i,1)-xmax)<tol && pos(i,2)>=ymin) ...
                || (abs(pos(i,2)-ymin)<tol && (pos(i,1) >= xmin && pos(i,1) <= xmax)) ...
                || (abs(pos(i,2)-ymax)<tol && (pos(i,1) >= xmin && pos(i,1) <= xmax))
            iatoms = iatoms + 1;
            interface_atom_list(iatoms,1) = i;
        end
    end
    
elseif strcmp(example,'Tensile_Test') || strcmp(example,'Tensile_Test_rotated')
    
    for i=1:length(pos(:,1))
        if (abs(pos(i,1)-xmin)<tol && pos(i,2)>=ymin) || (abs(pos(i,1)-xmax)<tol && pos(i,2)>=ymin)
            iatoms = iatoms + 1;
            interface_atom_list(iatoms,1) = i;
        end
    end
    
    
else
    error('Example not found')
end

interface_atom_list = unique(interface_atom_list);

end