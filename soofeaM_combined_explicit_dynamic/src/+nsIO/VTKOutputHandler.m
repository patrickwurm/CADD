classdef VTKOutputHandler < nsIO.OutputHandler
    %VTKOUTPUTHANDLER writes the results to a .vtk file
    
    properties ( SetAccess = protected )
    end
    
    methods
        function self = VTKOutputHandler( file_name )
            self@nsIO.OutputHandler( file_name );
        end
    end
    
    methods ( Static, Access = protected )
        function writeResults( model, file_id )            
            %Write file header
            fprintf(file_id,  '# vtk DataFile Version 3.0\n');
            fprintf(file_id,  'VTK Output\n');
            fprintf(file_id,  'ASCII\n');
            fprintf(file_id,  'DATASET UNSTRUCTURED_GRID\n');
            
            % Maybe put the following lines to a more appropriate
            % position...
            number_of_nodes = length(model.node_dict);
            fprintf(file_id,  ['\nPOINTS ',num2str(number_of_nodes),' float\n']);
            
            %Determine element type
            if isa(model.element_dict(1).type.shape,'nsModel.nsShape.QuadShape')
                element_type = 'quad';
            elseif isa(model.element_dict(1).type.shape,'nsModel.nsShape.HexShape')
                element_type = 'hex';
            elseif isa(model.element_dict(1).type.shape,'nsModel.nsShape.TriangleShape')
                element_type = 'triangle';
            elseif isa(model.element_dict(1).type.shape,'nsModel.nsShape.TetraShape')
                element_type = 'tetra';
            else
                model.element_dict(1).type.shape
                error(['VTK Output for element type "',class(model.element_dict(1).type.shape),'" not implemented.'])
            end
            
            switch element_type
                case 'quad'
                    %Write nodal coordinates
                    for node_counter = 1:number_of_nodes
                        out = [model.node_dict(node_counter).material_coordinates, 0]; %additional zero: there must be three coordinates!
                        fprintf(file_id, [num2str(out), '\n']);
                    end
                    %Write element information
                    number_of_elements = numel(model.element_dict);
                    fprintf(file_id,  ['\nCELLS ',num2str(number_of_elements),' ',num2str(5*number_of_elements),'\n']);
                    for element = model.element_dict
                        fprintf(file_id, ['4 ',num2str(element.node_number_list-1),'\n']);
                        %Node numbering for vtk starts at 0,
                        %that's why there is a '-1' in the line above.
                    end
                    fprintf(file_id, ['\nCELL_TYPES ',num2str(number_of_elements),'\n']);
                    for i=1:number_of_elements
                        fprintf(file_id,  '9\n');
                    end
                    %Write displacement vector
                    fprintf(file_id, ['\nPOINT_DATA ',num2str(number_of_nodes),'\n']);
                    fprintf(file_id,  'VECTORS displacement double\n');
                    for node_counter = 1:number_of_nodes
                        out = [model.node_dict(node_counter).dof.displacement, 0];%additional zero: there must be three displacements!
                        fprintf(file_id, [num2str(out), '\n']);
                    end
                    
                case 'hex'
                    %Write nodal coordinates
                    for node_counter = 1:number_of_nodes
                        out = [model.node_dict(node_counter).material_coordinates];
                        fprintf(file_id, [num2str(out), '\n']);
                    end
                    %Write element information
                    number_of_elements = numel(model.element_dict);
                    fprintf(file_id,  ['\nCELLS ' ,num2str(number_of_elements),' ',num2str(9*number_of_elements),'\n']);
                    for element = model.element_dict
                        fprintf(file_id, ['8 ',num2str(element.node_number_list-1), '\n']);
                        %Node numbering for vtk starts at 0,
                        %that's why there is a '-1' in the line above.
                    end
                    fprintf(file_id,  ['\nCELL_TYPES ',num2str(number_of_elements),'\n']);
                    for i=1:number_of_elements
                        fprintf(file_id,  '12\n');
                    end
                    %Write displacement vector
                    fprintf(file_id,  ['\nPOINT_DATA ',num2str(number_of_nodes), '\n']);
                    fprintf(file_id,  'VECTORS displacement double\n');
                    for node_counter = 1:number_of_nodes
                        out = model.node_dict(node_counter).dof.displacement;
                        fprintf(file_id, [num2str(out), '\n']);
                    end
                    
                case 'triangle'
                    %Write nodal coordinates
                    for node_counter = 1:number_of_nodes
                        out = [model.node_dict(node_counter).material_coordinates, 0]; %additional zero: there must be three coordinates!
                        fprintf(file_id, [num2str(out), '\n']);
                    end
                    %Write element information
                    number_of_elements = numel(model.element_dict);
                    fprintf(file_id,  ['\nCELLS ',num2str(number_of_elements),' ',num2str(4*number_of_elements),'\n']);
                    for element = model.element_dict
                        fprintf(file_id, ['3 ', num2str(element.node_number_list-1), '\n']);
                        %Node numbering for vtk starts at 0,
                        %that's why there is a '-1' in the line above.
                    end
                    fprintf(file_id,  ['\nCELL_TYPES ',num2str(number_of_elements), '\n']);
                    for i=1:number_of_elements
                        fprintf(file_id,  '5\n');
                    end
                    %Write displacement vector
                    fprintf(file_id,  ['\nPOINT_DATA ',num2str(number_of_nodes), '\n']);
                    fprintf(file_id,  'VECTORS displacement double\n');
                    for node_counter = 1:number_of_nodes
                        out = [model.node_dict(node_counter).dof.displacement, 0];%additional zero: there must be three displacements!
                        fprintf(file_id, [num2str(out), '\n']);
                    end
                case 'tetra'                
                    %Write nodal coordinates
                    for node_counter = 1:number_of_nodes
                        out = [model.node_dict(node_counter).material_coordinates];
                        fprintf(file_id, [num2str(out), '\n']);
                    end
                    %Write element information
                    number_of_elements = numel(model.element_dict);
                    fprintf(file_id,  ['\nCELLS ',num2str(number_of_elements),' ',num2str(5*number_of_elements),'\n']);
                    for element = model.element_dict
                        fprintf(file_id, ['4 ', num2str(element.node_number_list-1), '\n']);
                        %Node numbering for vtk starts at 0,
                        %that's why there is a '-1' in the line above.
                    end
                    fprintf(file_id,  ['\nCELL_TYPES ',num2str(number_of_elements), '\n']);
                    for i=1:number_of_elements
                        fprintf(file_id,  '10\n');
                    end
                    %Write displacement vector
                    fprintf(file_id,  ['\nPOINT_DATA ',num2str(number_of_nodes), '\n']);
                    fprintf(file_id,  'VECTORS displacement double\n');
                    for node_counter = 1:number_of_nodes
                        out = [model.node_dict(node_counter).dof.displacement];
                        fprintf(file_id, [num2str(out), '\n']);
                    end
            end
            
            % Output von mises stress
            fprintf(file_id,  '\nSCALARS stress_von_mises double 1');
            fprintf(file_id,  '\nLOOKUP_TABLE default\n');
            for node_counter = 1:number_of_nodes
                out = model.node_dict(node_counter).stress_von_mises;
                fprintf(file_id, [num2str(out), ' \n']);
            end
        end        
    end
end