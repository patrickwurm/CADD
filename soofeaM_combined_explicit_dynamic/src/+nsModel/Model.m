classdef Model < handle
    %MODEL The model of the finite element problem
    properties ( SetAccess = private )
        dimension
        time_step
        bc_handler = [];
        
        node_dict = nsModel.Node.empty
        element_dict = nsModel.Element.empty
        edge_dict = nsModel.Edge.empty
        face_dict = nsModel.Face.empty
        boundary_dict = nsModel.Boundary.empty
        dislocation_dict = nsModel.ContinuumDislocation.empty
        ghostdislocation_dict = nsModel.GhostDislocation.empty
        
        type_dict = {}
        material_dict = {}
        
        time_bar = [];
        last_finished_time_step
        example
        
        time %current atomistic time step
        a %atomic distance
    end
    
    methods
        function self = Model(dim)
            self.dimension = dim;
            self.last_finished_time_step = 0;
        end
        
        function setLastFinishedTimeStep(self, last_finished_time_step)
            self.last_finished_time_step = last_finished_time_step;
        end
        
        function addTimeBar(self, time_bar)
            self.time_bar = time_bar;
        end
        
        function setBCHandler( self, bc_handler )
            self.bc_handler = bc_handler;
            self.bc_handler.setModel(self);        
        end
        
        function setTimeStep( self, time_step )
            self.time_step = time_step;
        end
        
        function addNode(self, new_node)
            self.node_dict(new_node.number) = new_node;
        end
        
        function n_nodes = getNumberOfNodes(self)
            n_nodes = length(self.node_dict);
        end
        
        function addElement(self, new_element)
            node_list = self.resolveNodes( new_element.node_number_list );
            new_element.setNodeList(node_list);
            
            self.element_dict(new_element.number) = new_element;

            % We set the element type to be 1!
            new_element.setType( self.type_dict{1} );
            
            % We set the material type to be 1!
            new_element.setMaterial( self.material_dict{1} );
        end
        
        function addEdge( self, new_edge )
            node_list = self.resolveNodes( new_edge.node_number_list );
            new_edge.setNodeList(node_list);
            
            self.edge_dict(new_edge.number) = new_edge;
            
            % We set the edge type to be 2!
            new_edge.setType( self.type_dict{2} );
        end
        
        function addFace( self, new_face )
            node_list = self.resolveNodes( new_face.node_number_list );
            new_face.setNodeList(node_list);
            
            self.face_dict(new_face.number) = new_face;
                        
            % We set the face type to be 2!
            new_face.setType( self.type_dict{2} );
        end
        
        function addContinuumDislocation( self, new_dislocation )
            self.dislocation_dict(new_dislocation.number) = new_dislocation;
        end
        
        function addGhostDislocation( self, new_dislocation )
            self.ghostdislocation_dict(new_dislocation.number) = new_dislocation;
        end
        
        function addType( self, new_type )
            self.type_dict{new_type.number} = new_type;
        end
              
        function addMaterial( self, new_material )
            self.material_dict{new_material.number} = new_material;
        end
        
        function number_of_unknowns = getNumberOfUnknowns(self)
           % Works only for a one-field-formulation
           number_of_unknowns = self.dimension * numel(self.node_dict);
        end
               
        function appendComponentToBoundary( self, boundary_number, component_number, component_type )
            
            % The following condition is fulfilled for the first boundary.
            % The second one checks wether the current boundary number
            % already exists in the boundary dict.
            if isempty(self.boundary_dict) || ~ismember(boundary_number ,[self.boundary_dict.number])
                self.boundary_dict(boundary_number) = nsModel.Boundary(boundary_number);
            end
            
            % add component to the boundary
            boundary = findobj(self.boundary_dict, 'number', boundary_number);
            if strcmp( component_type, 'edge' )
                new_component = self.edge_dict(component_number);
                              
            elseif strcmp( component_type, 'face' ) 
                new_component = self.face_dict(component_number);
            end
            
            boundary.addComponent( new_component );
        end
        
        function setA(self, a)
            self.a = a;
        end
        
        function setTime(self, time)
            self.time = time;
        end
    end
        
    methods ( Access = private )
        function node_list = resolveNodes( self, node_number_list )
            % This method takes a list of node numbers and returns a list
            % of references in the same order to the corresponding node
            % objects.
            
            % node_number_list: The node numbers which should be resolved
            % https://de.mathworks.com/help/matlab/ref/empty.html
            len = length(node_number_list);
            node_list = nsModel.Node.empty(0,len);
            
            for i = 1:len
                node_list(i) = self.node_dict( node_number_list(i) );
            end
        end
    end
end