function findNatCoordsPad( fe_model )

%findNatCoordsPad uses the Locate Part of the "Generalized Particle
%Search-Locate Algorithm for Arbitrary Grids" by Allievi, Bermejo, 1997
%To find the natural coordinates of the Pad atoms.
%These coordinates are needed to find the pad atom displacements through
%interpolation.

disp('Computing Natural Coordinates of Pad Atoms ...')

global tol
conv_criteria = tol*1000;

for element = fe_model.element_dict
    a = element.node_list(1).material_coordinates;
    b = element.node_list(2).material_coordinates;
    c = element.node_list(3).material_coordinates;
    if ~isempty(element.pad_atoms)
        for ipa = 1:numel(element.pad_atoms) %loop over all pad atoms associated with the element
            xp = element.pad_atoms(ipa).coordinates(1);
            yp = element.pad_atoms(ipa).coordinates(2);
            r = 0;
            s = 0;
            CONV = false;
            while CONV == false
                delta = (b(1)-a(1))*(c(2)-a(2))-(c(1)-a(1))*(b(2)-a(2));
                X = element.getCoordinateArray('material');
                H = element.type.shape.getArray( [r, s] );
                Coords = transpose(X)*transpose(H);
                rold = r;
                sold = s;
                r = r + 1/delta * ( (c(2)-a(2)) * (xp - Coords(1)) +  (a(1)-c(1)) * (yp - Coords(2)) );
                s = s + 1/delta * ( (a(2)-b(2)) * (xp - Coords(1)) +  (b(1)-a(1)) * (yp - Coords(2)) );
                norm = sqrt( (r-rold)^2 + (s-sold)^2  );
                %plotelement(element, refpos, [r,s], element.pad_atoms(ipa))
                if norm < conv_criteria
                    CONV = true;
                    %throw an error if the natural coordinates are outside
                    %of the element (something must be wrong with the pad atom/element association)
                    if r <=-10^-12 || r >= 1+10^-12 || s <= -10^-12 || s >= 1+10^-12
                        error(['Pad atom ',num2str(element.pad_atoms(ipa).number),' is associated wrong.'])
                    end
                    element.pad_atoms(ipa).setNaturalCoords([r,s]); 
                end
            end
        end
    else %if there are no pad atoms associated with the element, do nothing.
    end
end

end

