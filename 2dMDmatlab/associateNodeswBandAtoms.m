function [ band_atom_ass_nodes ] = associateNodeswBandAtoms(band_atoms, pos, fe_model, xbandl, xbandr, yband, interface_node_list)

band_atom_ass_nodes = zeros(size(band_atoms));

%seperate the band atoms into two groups, the horizontal ones and the
%vertical ones
ba_vert = [];
ba_hor = [];

for i = 1:length(band_atoms)
    if pos(band_atoms(i),2) >= yband + 10^-15
        ba_vert(end+1) = band_atoms(i);
    else
        ba_hor(end+1) = band_atoms(i);
    end
end

%associate vertical nodes
for i = ba_vert
    associated = false;
    for j = interface_node_list'
        if abs(fe_model.node_dict(j).material_coordinates(2)  - pos(i,2)) <= 10^-15 && sum(band_atom_ass_nodes==j)==0 %this is the check if not already included
            band_atom_ass_nodes(find(band_atoms == i)) = j;
            associated = true;
            break
        end
    end
    if associated == false
        error(['Vertical node', num2str(i),' is not associated with any interface node.'])
    end
end

%associate horizontal nodes
for i = ba_hor
    for j = interface_node_list'
        if abs(fe_model.node_dict(j).material_coordinates(1)  - pos(i,1)) <= 10^-15 && sum(band_atom_ass_nodes==j)==0 %this is the check if not already included
            band_atom_ass_nodes(find(band_atoms == i)) = j;
        end
    end
    if associated == false
        error(['Horizontal node', num2str(i),' is not associated with any interface node.'])
    end
end

%UNIT TEST
%TEST_associateNodeswBandAtoms(band_atom_ass_nodes, band_atoms, pos, fe_model)


end

