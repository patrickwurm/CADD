function dump = getPosDump(example_name, mstep, refpos, pos, fe_model, fe_analysis, istep, dump)

    if strcmp(example_name, 'Example2')  || strcmp(example_name, 'Example2_rotated')
        dump.group1a_positions(istep+1) = mean(pos(dump.group1a_atoms,1));
        dump.group1b_positions(istep+1) = mean(pos(dump.group1b_atoms,1));
        dump.group1c_positions(istep+1) = mean(pos(dump.group1c_atoms,1));
        
        dump.group2a_positions(istep+1) = mean(pos(dump.group2a_atoms,1));
        dump.group2b_positions(istep+1) = mean(pos(dump.group2b_atoms,1));
        dump.group2c_positions(istep+1) = mean(pos(dump.group2c_atoms,1));
        
        dump.group3a_positions(istep+1) = mean(pos(dump.group3a_atoms,1));
        dump.group3b_positions(istep+1) = mean(pos(dump.group3b_atoms,1));
        dump.group3c_positions(istep+1) = mean(pos(dump.group3c_atoms,1));
        
        for i = 1:length(dump.node1_atoms)
            node = fe_model.node_dict(1,dump.node1_atoms(i));
            global_index = (node.number-1) * fe_model.dimension + 1;
            dump.node1_positions(istep+1) = dump.node1_positions(istep+1) + node.material_coordinates(1) + fe_analysis.u(global_index);
        end
        for i = 1:length(dump.node2_atoms)
            node = fe_model.node_dict(1,dump.node2_atoms(i));
            global_index = (node.number-1) * fe_model.dimension + 1;
            dump.node2_positions(istep+1) = dump.node2_positions(istep+1) + node.material_coordinates(1) + fe_analysis.u(global_index);
        end
        for i = 1:length(dump.node3_atoms)
            node = fe_model.node_dict(1,dump.node3_atoms(i));
            global_index = (node.number-1) * fe_model.dimension + 1;
            dump.node3_positions(istep+1) = dump.node3_positions(istep+1) + node.material_coordinates(1) + fe_analysis.u(global_index);
        end
        dump.node1_positions(istep+1) = dump.node1_positions(istep+1)/length(dump.node1_atoms);
        dump.node2_positions(istep+1) = dump.node2_positions(istep+1)/length(dump.node2_atoms);
        dump.node3_positions(istep+1) = dump.node3_positions(istep+1)/length(dump.node3_atoms);
    else
        %do nothing
    end

end

