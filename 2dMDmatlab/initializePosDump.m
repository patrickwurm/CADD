function dump = initializePosDump(example_name, mstep, refpos, fe_model)

dump = [];

if strcmp(example_name, 'Example2') || strcmp(example_name, 'Example2_rotated')
    dump.group1a_positions = zeros(1,mstep+1);
    dump.group1b_positions = zeros(1,mstep+1);
    dump.group1c_positions = zeros(1,mstep+1);
    dump.group2a_positions = zeros(1,mstep+1);
    dump.group2b_positions = zeros(1,mstep+1);
    dump.group2c_positions = zeros(1,mstep+1);
    dump.group3a_positions = zeros(1,mstep+1);
    dump.group3b_positions = zeros(1,mstep+1);
    dump.group3c_positions = zeros(1,mstep+1);
    dump.node1_positions = zeros(1,mstep+1);
    dump.node2_positions = zeros(1,mstep+1);
    dump.node3_positions = zeros(1,mstep+1);
    
    dump.group1a_atoms = [];
    dump.group1b_atoms = [];
    dump.group1c_atoms = [];
    dump.group2a_atoms = [];
    dump.group2b_atoms = [];
    dump.group2c_atoms = [];
    dump.group3a_atoms = [];
    dump.group3b_atoms = [];
    dump.group3c_atoms = [];
    dump.node1_atoms = [];
    dump.node2_atoms = [];
    dump.node3_atoms = [];
    
    a=fe_model.a;
    postol = 0.01*a;
    mxatm = length(refpos(:,1));
    for i = 1:mxatm
        if refpos(i,1) >= 10*a-postol && refpos(i,1) <= 15*a+postol
            dump.group1a_atoms(end+1) = i;
        end
        if refpos(i,1) >= 20*a-postol && refpos(i,1) <= 25*a+postol
            dump.group1b_atoms(end+1) = i;
        end
        if refpos(i,1) >= 30*a-postol && refpos(i,1) <= 35*a+postol
            dump.group1c_atoms(end+1) = i;
        end
        if refpos(i,1) >= 90*a-postol && refpos(i,1) <= 95*a+postol
            dump.group2a_atoms(end+1) = i;
        end
        if refpos(i,1) >= 100*a-postol && refpos(i,1) <= 105*a+postol
            dump.group2b_atoms(end+1) = i;
        end
        if refpos(i,1) >= 110*a-postol && refpos(i,1) <= 115*a+postol
            dump.group2c_atoms(end+1) = i;
        end
        if refpos(i,1) >= 165*a-postol && refpos(i,1) <= 170*a+postol
            dump.group3a_atoms(end+1) = i;
        end
        if refpos(i,1) >= 175*a-postol && refpos(i,1) <= 180*a+postol
            dump.group3b_atoms(end+1) = i;
        end
        if refpos(i,1) >= 185*a-postol && refpos(i,1) <= 190*a+postol
            dump.group3c_atoms(end+1) = i;
        end
    end
    
    for node = fe_model.node_dict
        noderefpos = node.material_coordinates(1);
        if abs(noderefpos-(-5*a))<=1e-15
            dump.node1_atoms(end+1) = node.number; %not really atoms but nodes
        end
        if abs(noderefpos-(-95.3125*a))<=a
            dump.node2_atoms(end+1) = node.number; %not really atoms but nodes
        end
        if abs(noderefpos-(-171.25*a))<=a
            dump.node3_atoms(end+1) = node.number; %not really atoms but nodes
        end
    end
    
    dump.group1a_positions(1) = mean(refpos(dump.group1a_atoms,1));
    dump.group1b_positions(1) = mean(refpos(dump.group1b_atoms,1));
    dump.group1c_positions(1) = mean(refpos(dump.group1c_atoms,1));
    
    dump.group2a_positions(1) = mean(refpos(dump.group2a_atoms,1));
    dump.group2b_positions(1) = mean(refpos(dump.group2b_atoms,1));
    dump.group2c_positions(1) = mean(refpos(dump.group2c_atoms,1));
    
    dump.group3a_positions(1) = mean(refpos(dump.group3a_atoms,1));
    dump.group3b_positions(1) = mean(refpos(dump.group3b_atoms,1));
    dump.group3c_positions(1) = mean(refpos(dump.group3c_atoms,1));
    
    for i = 1:length(dump.node1_atoms)
        dump.node1_positions(1) = dump.node1_positions(1) + fe_model.node_dict(1,dump.node1_atoms(i)).material_coordinates(1);
    end
    for i = 1:length(dump.node2_atoms)
        dump.node2_positions(1) = dump.node2_positions(1) + fe_model.node_dict(1,dump.node2_atoms(i)).material_coordinates(1);
    end
    for i = 1:length(dump.node3_atoms)
        dump.node3_positions(1) = dump.node3_positions(1) + fe_model.node_dict(1,dump.node3_atoms(i)).material_coordinates(1);
    end
    dump.node1_positions(1) = dump.node1_positions(1)/length(dump.node1_atoms);
    dump.node2_positions(1) = dump.node2_positions(1)/length(dump.node2_atoms);
    dump.node3_positions(1) = dump.node3_positions(1)/length(dump.node3_atoms);
    
end

end

