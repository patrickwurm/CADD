classdef GmshInputHandler < nsIO.InputHandler
    %GMSHINPUTHANDLER reads Gmsh generated '.msh' files
    
    properties ( SetAccess = protected )
        element_types
        meshes = {}
        number_of_nodes = 0
        number_of_elements = 0
    end
    
    methods
        function self = GmshInputHandler( file_name, dimension )
            self@nsIO.InputHandler( file_name, dimension );
            
            gmsh_element_numbers = { '1', '2', '3', '4', '5', '15' };
            gmsh_element_types = { '2-node_line', '3-node_triangle', '4-node_quadrangle', '4-node_tetrahedron', '8-node_hexahedron', '1-node_point' };
            self.element_types = containers.Map( gmsh_element_numbers, gmsh_element_types );
        end
        
        function readMesh ( self, inp_model, file_id )
            disp('-> Parsing Gmsh input file...');
            state = 'undefined';
            boundary_ids = [];
            
            % These variables are updated in the regarding sections
            number_of_points = 0;
            number_of_edges = 0;
            number_of_faces = 0;
            
            current_line = 'start';
            while ischar( current_line )
                current_line = fgetl( file_id );
                % Parsing is controlled via the state-variable.
                if( strcmp(state,'undefined') )
                    if( strcmp(current_line,'$PhysicalNames') )
                        state = 'physical';
                        fprintf('--> Reading physical names...\n');
                        continue;
                    end
                    if( strcmp(current_line,'$Nodes') )
                        state = 'nodes';
                        fprintf('--> Reading nodes...\n');
                        continue;
                    end
                    if( strcmp(current_line,'$Elements') )
                        state = 'elements_container';
                        fprintf('--> Reading element container...\n');
                        continue;
                    end
                end
                
                % Reading physical groups. Note that in general all
                % physical entities appear here (e.g. physical lines). Yet
                % we only consider the case of physical entities
                % representing a mesh.
                if( strcmp(state,'physical') )
                    
                    % When end of PhysicalNames-container is reached, set
                    % back to 'undefined' and go to next line.
                    if( strcmp(current_line,'$EndPhysicalNames' ) )
                        state = 'undefined';
                        continue
                    end
                    
                    % Splits the string containing the current line at the
                    % spaces into a cell array.
                    line_data = strsplit( current_line );
                    
                    % We don't care about the number of named physical
                    % entities.
                    if( numel(line_data) == 1 )
                        continue
                    end
                    
                    % The name of the physical entity is in the third
                    % column. Remove the double quotes.
                    mesh_name = strrep(line_data{3},'"','');
                    
                    % check wether incoming mesh is new
                    mesh_is_new = true;
                    for i = 1 : numel(self.meshes)
                        if strcmp( self.meshes{i}, mesh_name )
                            mesh_is_new = false;
                        end
                    end
                    
                    if mesh_is_new
                        self.meshes = [self.meshes mesh_name];
                    end
                end
                
                % Reading nodes
                if( strcmp(state,'nodes') )
                    
                    % When end of nodes container is reached, set state
                    % back to 'undefined' and go to next line.
                    if( strcmp(current_line,'$EndNodes' ) )
                        state = 'undefined';
                        continue
                    end
                    
                    % Splits the string containing the current line at the
                    % spaces into a cell array.
                    line_data = strsplit( current_line );
                    
                    % The first line in the node container holds the number
                    % of nodes in the mesh.
                    if( numel(line_data) == 1 )
                        self.number_of_nodes = line_data{1};
                        continue
                    end
                    
                    switch self.dimension
                        case 1
                            new_node = nsModel.Node( str2double(line_data{1}), [ str2double(line_data{2}) ] );
                            inp_model.addNode( new_node );
                            continue
                        case 2
                            new_node = nsModel.Node( str2double(line_data{1}), [ str2double(line_data{2}) str2double(line_data{3}) ] );
                            inp_model.addNode( new_node );
                            continue
                        case 3
                            new_node = nsModel.Node( str2double(line_data{1}), [ str2double(line_data{2}) str2double(line_data{3}) str2double(line_data{4}) ] );
                            inp_model.addNode( new_node );
                            continue
                        otherwise
                            % Maybe the dimension should be checked in the
                            % constructor
                            assert(false,'Dimension must be 1,2 or 3.')
                    end
                end
                
                % Reading element container
                if( strcmp(state,'elements_container') || strcmp(state,'elements') || strcmp(state,'faces')|| strcmp(state,'edges') || strcmp(state,'points') )
                    
                    % When end of elements container is reached, set state
                    % back to 'undefined' and go to next line.
                    if( strcmp(current_line,'$EndElements' ) )
                        state = 'undefined';
                        continue
                    end
                    
                    % Splits the string containing the current line at the
                    % spaces into a cell array.
                    line_data = strsplit( current_line );
                    
                    % The first line in the elements container holds the number
                    % of elements in the mesh.
                    if( numel(line_data) == 1 )
                        self.number_of_elements = line_data{1};
                        continue
                    end
                    
                    % detect if we are at the points, edges, faces or
                    % elements.
                    state = self.elementContainerState( line_data , state );
                    
                    % At the moment, we are not interested in points. they
                    % would be faces in a 1D-computation
                    if( strcmp(state,'points') )
                        number_of_points = number_of_points +1;
                        continue
                    end
                    
                    if( strcmp(state,'edges') )
                        
                        % Check number of tags in the msh-file
                        number_of_tags = str2double(line_data{3});
                        assert( number_of_tags==2, 'We need the id of the line were the edge lies on!' )
                        
                        % Start- and endpoint of the node numbers in the
                        % msh-file.
                        start_nodes = 4 + number_of_tags;
                        end_nodes = numel(line_data);
                        number_of_edge_nodes = end_nodes - start_nodes +1;
                        
                        node_number_list = zeros( 1, number_of_edge_nodes );
                        for i = 1 : number_of_edge_nodes
                            node_number_list(i) = str2double( line_data{ i+start_nodes-1 } );
                        end
                        
                        % add edge to the model. in the msh-file,
                        % the edges start after the
                        % points; hence the substraction.
                        edge_number = str2double(line_data{1}) - number_of_points;
                        inp_model.addEdge( nsModel.Edge( edge_number, node_number_list ) );
                        number_of_edges = number_of_edges + 1;
                        
                        % get boundary_id of current edge. This is the
                        % fourth entry in the msh-file.
                        new_boundary_id = str2double( line_data{4} );
                        
                        % add current boundary id to the boundary_ids-array
                        % if it doesn't exist yet.
                        if not( any(boundary_ids == new_boundary_id) )
                            boundary_ids = [boundary_ids new_boundary_id];
                        end
                        
                        if self.dimension == 2
                            inp_model.appendComponentToBoundary( new_boundary_id, edge_number, 'edge' );
                        end
                    end
                    
                    if( strcmp(state,'faces') )
                        
                        % Check number of tags in the msh-file
                        number_of_tags = str2double(line_data{3});
                        assert( number_of_tags==2, 'We need the id of the surface the face belongs to!' )
                        
                        % Start- and endpoint of the node numbers in the
                        % msh-file.
                        start_nodes = 4 + number_of_tags;
                        end_nodes = numel(line_data);
                        number_of_face_nodes = end_nodes - start_nodes +1;
                        
                        node_number_list = zeros( 1, number_of_face_nodes );
                        for i = 1 : number_of_face_nodes
                            node_number_list(i) = str2double( line_data{ i+start_nodes-1 } );
                        end
                        
                        % add face to the model. in the msh-file,
                        % the faces start after the
                        % edges; hence the substraction.
                        face_number = str2double(line_data{1}) - ( number_of_edges + number_of_points );
                        inp_model.addFace( nsModel.Face( face_number, node_number_list ) );
                        number_of_faces = number_of_faces + 1;
                        
                        % get boundary_id of current face. This is the
                        % fourth entry in the msh-file.
                        new_boundary_id = str2double( line_data{4} );
                        
                        % add current boundary id to the boundary_ids-array
                        % if it doesn't exist yet.
                        if not( any(boundary_ids == new_boundary_id) )
                            boundary_ids = [boundary_ids new_boundary_id];
                        end
                        
                        if self.dimension == 3
                            inp_model.appendComponentToBoundary( new_boundary_id, face_number, 'face' );
                        end
                    end
                    
                    if( strcmp(state,'elements') )
                        % Check number of tags in the msh-file
                        number_of_tags = str2double(line_data{3});
                        assert( number_of_tags==2, 'We need the id of the surface the face belongs to!' )
                        
                        % Start- and endpoint of the node numbers in the
                        % msh-file.
                        start_nodes = 4 + number_of_tags;
                        end_nodes = numel(line_data);
                        number_of_element_nodes = end_nodes - start_nodes +1;
                        
                        node_number_list = zeros( 1, number_of_element_nodes );
                        for i = 1 : number_of_element_nodes
                            node_number_list(i) = str2double( line_data{ i+start_nodes-1 } );
                        end
                        
                        % add element to the model. in the msh-file,
                        % the elements start after the
                        % faces; hence the substraction.
                        element_number = str2double(line_data{1}) - ( number_of_faces + number_of_edges + number_of_points ) ;
                        new_element = nsModel.Element( element_number, node_number_list );
                        inp_model.addElement( new_element );
                    end
                end
            end
        end
    end
    
    methods ( Access = private )
        function state = elementContainerState( self, line_data , old_state )
            
            % line_data{2}: the second column in the msh-file contains the
            % element type.
            if strcmp( self.element_types(line_data{2}), '1-node_point' )
                state = 'points';
            elseif strcmp( self.element_types(line_data{2}), '3-node_triangle' )
                if self.dimension == 2
                    state = 'elements';
                else
                    state = 'faces';
                end
            elseif strcmp( self.element_types(line_data{2}), '2-node_line' )
                state = 'edges';
            elseif strcmp( self.element_types(line_data{2}), '4-node_quadrangle' )
                if self.dimension == 2
                    state = 'elements';
                else
                    state = 'faces';
                end
            elseif strcmp( self.element_types(line_data{2}), '4-node_tetrahedron' )
                state = 'elements';
            elseif strcmp( self.element_types(line_data{2}), '8-node_hexahedron' )
                state = 'elements';
            else
                assert(false, 'We do not know element type : ')
            end
            
            % Output if the state changes
            if not( strcmp(state,old_state) )
                fprintf(['---> We are on new state ' state '\n']);
            end
        end
    end
end