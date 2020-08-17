classdef VTKCADDOutputHandler < nsIO.OutputHandler
    %VTKOUTPUTHANDLER writes the results to a .vtk file
    
    properties ( SetAccess = protected )
    end
    
    methods
        function self = VTKCADDOutputHandler( file_name, folder_name, dimension, force_overwrite )
            self@nsIO.OutputHandler( file_name, folder_name, dimension, force_overwrite );
        end
    end
    
    methods ( Static )
        function writeResults( model, file_id, atomrefpos, atompos, mxatm, mxpatm, step, energies, atom_velocities, nodal_velocities )
             number_of_atoms = length(atomrefpos(:,1));
            %Determine element type
            if isa(model.element_dict(1).type.shape,'nsModel.nsShape.TriangleShape')
                element_type = 'triangle';
            else
                model.element_dict(1).type.shape
                error(['VTK Output for element type "',class(model.element_dict(1).type.shape),'" not implemented.'])
            end
            
            switch element_type
                case 'triangle'
                    %Write displacement vector
                    fprintf(file_id,  strjoin({'POINT_DATA',num2str(model.getNumberOfNodes()+number_of_atoms)}));
                    fprintf(file_id,  '\n');
                    fprintf(file_id,  'VECTORS displacement double\n');
                    for node_counter = 1:model.getNumberOfNodes()
                        out = [model.node_dict(node_counter).dof.displacement, 0];%additional zero: there must be three displacements!
                        fprintf(file_id, '%.15e ', out);
                        fprintf(file_id, '\n');
                    end
                    for ia = 1:number_of_atoms
                        out = [atompos(ia,:)-atomrefpos(ia,:), 0];%additional zero: there must be three displacements!
                        fprintf(file_id, '%.15e ', out);
                        fprintf(file_id, '\n');
                    end
                    
                    %Write velocity vector
                    fprintf(file_id,  ' \n');
                    fprintf(file_id,  'VECTORS velocity double\n');
                    for node_counter = 1:model.getNumberOfNodes()
                        out = zeros(1,3);
                        for coord_index = 1:2
                            global_index = (node_counter-1) * 2 + coord_index;
                            out(coord_index) = nodal_velocities(global_index);
                        end
                        fprintf(file_id, '%.15e ', out);
                        fprintf(file_id, '\n');
                    end
                    atom_velocities = vertcat(atom_velocities, zeros(mxpatm,2));
                    for ia = 1:number_of_atoms
                        out = [atom_velocities(ia,:), 0];%additional zero: there must be three displacements!
                        fprintf(file_id, '%.15e ', out);
                        fprintf(file_id, '\n');
                    end
            end
            
            fprintf(file_id,  ' \n');
            fprintf(file_id,  'SCALARS stress_vonmises double 1\n');
            fprintf(file_id,  'LOOKUP_TABLE default\n');
            for node_counter = 1:model.getNumberOfNodes()
                if model.last_finished_time_step == 0 %only write undeformed configuration to file
                    out = 0;
                else
                    if isempty(model.node_dict(node_counter).stress_vonmises)
                        out = 0;
                    else
                        out = model.node_dict(node_counter).stress_vonmises;
                    end
                end
                fprintf(file_id, '%.15e ', out);
                fprintf(file_id, '\n');
            end
            for node_counter = 1:number_of_atoms
                out = 0;
                fprintf(file_id, '%.15e ', out);
                fprintf(file_id, '\n');
            end
            fprintf(file_id,  ' \n');
            fprintf(file_id,  ' \n');
            fprintf(file_id,  'SCALARS types double 1\n');
            fprintf(file_id,  'LOOKUP_TABLE default\n');
            for node_counter = 1:model.getNumberOfNodes()
                fprintf(file_id, '1\n');
            end
            for node_counter = 1:number_of_atoms
                 fprintf(file_id, '-1\n');
            end
            fprintf(file_id, ' \n');
            fprintf(file_id,  ' \n');
            fprintf(file_id,  'SCALARS energies double 1\n');
            fprintf(file_id,  'LOOKUP_TABLE default\n');
            for node_counter = 1:model.getNumberOfNodes()
                out = 0;
                fprintf(file_id, '%.15e ', out);
                fprintf(file_id, '\n');
            end
            for node_counter = 1:number_of_atoms
                out = energies(node_counter);
                fprintf(file_id, '%.15e ', out);
                fprintf(file_id, '\n');
            end
            fprintf(file_id,  ' \n');

        end
        
        function writeResults_standard( model, file_id, atomrefpos, atompos, mxatm, mxpatm, step )
            number_of_atoms = length(atomrefpos(:,1));
            %Write file header
            fprintf(file_id,  '# vtk DataFile Version 3.0\n');
            fprintf(file_id,  ['mxatm ',num2str(mxatm), ' mxpatm ', num2str(mxpatm), ' step ', num2str(step),'\n']);
            fprintf(file_id,  'ASCII\n');
            fprintf(file_id,  'DATASET UNSTRUCTURED_GRID\n');
            fprintf(file_id,  [strjoin({'POINTS',num2str(model.getNumberOfNodes()+number_of_atoms),'float'}),'\n']);
            
            %Determine element type
            if isa(model.element_dict(1).type.shape,'nsModel.nsShape.TriangleShape')
                element_type = 'triangle';
            else
                model.element_dict(1).type.shape
                error(['VTK Output for element type "',class(model.element_dict(1).type.shape),'" not implemented.'])
            end
            
            switch element_type
                case 'triangle'
                    %Write nodal coordinates
                    for node_counter = 1:model.getNumberOfNodes()
                        out = [model.node_dict(node_counter).spatial_coordinates, 0]; %additional zero: there must be three coordinates!
                        fprintf(file_id, '%.15e ', out);
                        fprintf(file_id, '\n');
                    end
                    for ia = 1:number_of_atoms
                        out = [atomrefpos(ia,:), 0]; %additional zero: there must be three coordinates!
                        fprintf(file_id, '%.15e ', out);
                        fprintf(file_id, '\n');
                    end
                    
                    fprintf(file_id, '\n');
                    %Write element information
                    number_of_elements = numel(model.element_dict);
                    fprintf(file_id,  strjoin({'CELLS',num2str(number_of_elements+number_of_atoms),num2str(4*number_of_elements+2*number_of_atoms)}));
                    fprintf(file_id, '\n');
                    for element = model.element_dict
                        fprintf(file_id,  [num2str(3),' ',num2str(element.node_number_list-1),'\n']);
                        %Node numbering for vtk starts at 0,
                        %that's why there is a '-1' in the line above.
                    end
                    for ia = 1:number_of_atoms
                         fprintf(file_id, [num2str(1),' ',num2str(model.getNumberOfNodes()+ia-1),'\n']);
                    end
                    fprintf(file_id, ' \n');
                    fprintf(file_id,  strjoin({'CELL_TYPES',num2str(number_of_elements+number_of_atoms)}));
                    fprintf(file_id,  '\n');
                    for i=1:number_of_elements
                        fprintf(file_id,  '5\n');
                    end
                    for i=1:number_of_atoms
                        fprintf(file_id,  '1\n');
                    end
                    fprintf(file_id,  ' \n');
            end
        end
        
    end
end
