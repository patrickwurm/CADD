function interface_node_list = associateIntAtomswNodes(fe_model, refpos, interface_atom_list)

disp(' ')
disp('Associating Interface Atoms with Nodes ...')
interface_node_list = zeros(size(interface_atom_list));
global tol

for i = 1:numel(interface_atom_list)
    found_node = false;
    int_atom_nr = interface_atom_list(i);
    for node = fe_model.node_dict
        if abs((node.material_coordinates(1)-refpos(int_atom_nr,1)))<tol && abs((node.material_coordinates(2)-refpos(int_atom_nr,2)))<tol
            found_node = true;
            node.setInterfaceAtomNumber(int_atom_nr);
            interface_node_list(i) = node.number;
            break
        end
    end
    if found_node == false
        error(['Could not find node corresponding to coordinates: ',num2str(refpos(int_atom_nr,1)),' ',num2str(refpos(int_atom_nr,2))])
    end
    
end

end