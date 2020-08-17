function [ band_atom_ass_nodes ] = associateNodeswBandAtoms(example_name, band_atoms, pos, fe_model, xbandl, xbandr, ybandb, ybandt, interface_atom_list)

if strcmp(example_name, 'Radial_Pulse')
    band_atom_ass_nodes = zeros(size(band_atoms));
    
    %seperate the band atoms into two groups, the horizontal ones and the
    %vertical ones
    ba_vert = [];
    ba_hor = [];
    
    for i = 1:length(band_atoms)
        if pos(band_atoms(i),2) > ybandb + 10^-15 &&  pos(band_atoms(i),2) < ybandt - 10^-15
            ba_vert(end+1) = band_atoms(i);
        else
            ba_hor(end+1) = band_atoms(i);
        end
    end
    
    %associate vertical nodes
    for i = ba_vert
        associated = false;
        k = 0;
        for j = interface_atom_list'
            if abs(pos(j,2)  - pos(i,2)) <= 10^-15 %&& sum(band_atom_ass_nodes==j)==0 %this is the check if not already included
                k=k+1;
                two_atoms(k) = j;
                if k==2
                    if abs(pos(two_atoms(1),1)  - pos(i,1)) < abs(pos(two_atoms(2),1)  - pos(i,1))
                        band_atom_ass_nodes(find(band_atoms == i)) = two_atoms(1);
                        associated = true;
                        %                     scatter(pos(i,1),pos(i,2))
                        %                     hold on
                        %                     scatter(pos(two_atoms(1),1), pos(two_atoms(1),2))
                        break
                    else
                        band_atom_ass_nodes(find(band_atoms == i)) = two_atoms(2);
                        associated = true;
                        %                     scatter(pos(i,1),pos(i,2))
                        %                     hold on
                        %                     scatter(pos(two_atoms(2),1), pos(two_atoms(2),2))
                        break
                    end
                end
            end
        end
        if associated == false
            error(['Vertical node ', num2str(i),' is not associated with any interface node.'])
        end
    end
    
    %associate horizontal nodes
    for i = ba_hor
        associated = false;
        k=0;
        for j = interface_atom_list'
            if abs(pos(j,1)  - pos(i,1)) <= 10^-15 %&& sum(band_atom_ass_nodes==j)==0 %this is the check if not already included
                k=k+1;
                two_atoms(k) = j;
                if k==2
                    if abs(pos(two_atoms(1),2)  - pos(i,2)) < abs(pos(two_atoms(2),2)  - pos(i,2))
                        band_atom_ass_nodes(find(band_atoms == i)) = two_atoms(1);
                        associated = true;
                        %                     scatter(pos(i,1),pos(i,2))
                        %                     hold on
                        %                     scatter(pos(two_atoms(1),1), pos(two_atoms(1),2))
                        break
                    else
                        band_atom_ass_nodes(find(band_atoms == i)) = two_atoms(2);
                        associated = true;
                        %                     scatter(pos(i,1),pos(i,2))
                        %                     hold on
                        %                     scatter(pos(two_atoms(2),1), pos(two_atoms(2),2))
                        break
                    end
                end
            end
        end
        if associated == false
            error(['Horizontal node', num2str(i),' is not associated with any interface node.'])
        end
    end
elseif strcmp(example_name, 'Example2')
    
    band_atom_ass_nodes = zeros(size(band_atoms));

    %associate vertical nodes
    for i = band_atoms
        associated = false;
        for j = interface_atom_list'
            if abs(pos(j,2)  - pos(i,2)) <= 10^-15 %&& sum(band_atom_ass_nodes==j)==0 %this is the check if not already included
                        band_atom_ass_nodes(find(band_atoms == i)) = j;
                        associated = true;
                        break
            end
        end
        if associated == false
            error(['Vertical node ', num2str(i),' is not associated with any interface node.'])
        end
    end
elseif strcmp(example_name, 'Example2_rotated')
    
    band_atom_ass_nodes = zeros(size(band_atoms));

    %associate vertical nodes
    for i = band_atoms
        associated = false;
        for j = interface_atom_list'
            if abs(pos(j,2)  - pos(i,2)) <= 10^-15 %&& sum(band_atom_ass_nodes==j)==0 %this is the check if not already included
                        band_atom_ass_nodes(find(band_atoms == i)) = j;
                        associated = true;
                        break
            end
        end
        if associated == false
            error(['Vertical node ', num2str(i),' is not associated with any interface node.'])
        end
    end
    
elseif strcmp(example_name, 'Tensile_Test') || strcmp(example_name, 'Tensile_Test_rotated')
    band_atom_ass_nodes = zeros(size(band_atoms));
    
    %seperate the band atoms into two groups, the horizontal ones and the
    %vertical ones

    % figure(77)
    % scatter(pos(band_atoms,1),pos(band_atoms,2),'k')
    % hold on
    % scatter(pos(interface_atom_list,1),pos(interface_atom_list,2),'k')
    
    %associate vertical nodes
    for i = band_atoms
        associated = false;
        k = 0;
        for j = interface_atom_list'
            if abs(pos(j,2)  - pos(i,2)) <= 10^-15 %&& sum(band_atom_ass_nodes==j)==0 %this is the check if not already included
                k=k+1;
                two_atoms(k) = j;
                if k==2
                    if abs(pos(two_atoms(1),1)  - pos(i,1)) < abs(pos(two_atoms(2),1)  - pos(i,1))
                        band_atom_ass_nodes(find(band_atoms == i)) = two_atoms(1);
                        associated = true;
                        %                     scatter(pos(i,1),pos(i,2))
                        %                     hold on
                        %                     scatter(pos(two_atoms(1),1), pos(two_atoms(1),2))
                        break
                    else
                        band_atom_ass_nodes(find(band_atoms == i)) = two_atoms(2);
                        associated = true;
                        %                     scatter(pos(i,1),pos(i,2))
                        %                     hold on
                        %                     scatter(pos(two_atoms(2),1), pos(two_atoms(2),2))
                        break
                    end
                end
            end
        end
        if associated == false
            error(['Vertical node ', num2str(i),' is not associated with any interface node.'])
        end
    end

else
    error('Example not implemented.')
end

%UNIT TEST
%TEST_associateNodeswBandAtoms(band_atom_ass_nodes, band_atoms, pos, interface_atom_list)
end

