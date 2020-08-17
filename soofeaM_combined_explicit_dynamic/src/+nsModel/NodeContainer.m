classdef NodeContainer < nsModel.NumberedObject
    %NODECONTAINER Contains Nodes - superclass for Edge, Face and Element
    
    properties ( SetAccess = protected )
        node_list = nsModel.Node.empty
        node_number_list
        type
        int_points
    end
    
    methods
        function self = NodeContainer(number, node_number_list)
            self@nsModel.NumberedObject( number );
            
            % The list of node numbers.  This is redundant as we also store
            % 'node_list' as a property. But the problem is, that at the moment of
            % construction (during reading the input file) the node references do not
            % exist. Therefore we first provide the node number list and then resolve
            % this list to real node references in soofea.model.Model.resolveNodes().
            self.node_number_list = node_number_list;
            
            % An array containing the node references. It is important to
            % note that the nodes inside the 'node_list' are sorted
            % regarding their local node numbers.
            self.node_list = [];
            
            % Containing the model.Type of this
            % object
            self.type = [];
            
            % Containing all the soofea.model.RealIntegrationPoint of this
            % object. This list is initialized by self.setType as soon as a
            % type is assigned.
            self.int_points = [];
        end
               
        function setNodeList(self, node_list)
            self.node_list = node_list;
        end
        
        function number_of_nodes = getNumberOfNodes(self)
            number_of_nodes = length(self.node_list);
        end
        
        % liefert Koordinaten der Knoten des NodeContainers im Format:
        %
        %   [ x1 y1 z1;
        %     x2 y2 z2;
        %       ...
        %     xn yn zn ]
        %
        % 'n' ist dabei die Anzahl der Knoten des NodeContainers.
        function coordinate_array = getCoordinateArray( self, configuration )
            % In the linear case, the configuration is always 'material'.
            % Purtroppo, non ci sono dei parametri di default in matlab.
            % Ich glaube es ist besser, wenn man die Konfiguration angeben MUSS,
            % damit man in der nonlinear-implementation nichts übersieht.
            % dadurch steht zwar im linear case unnötigerweise 'material',
            % aber was solls...
%             if nargin == 1
%                 configuration = 'material';
%             end
            
            number_of_nodes = self.getNumberOfNodes();
            dimension = self.node_list(1).dimension;
            
            % coordinate_array is an ( number_of_nodes x dimension ) array.
            coordinate_array = zeros( number_of_nodes, dimension );
            
            if strcmp(configuration, 'material')
                for i = 1:number_of_nodes
                    coordinate_array( i, : ) = self.node_list(i).material_coordinates;
                end
            elseif strcmp(configuration, 'spatial')
                for i = 1:number_of_nodes
                    coordinate_array( i, : ) = self.node_list(i).spatial_coordinates;
                end
            else
                error('Configuration needs to be string ''material'' or ''spatial''!')
            end
        end
        
        % liefert Knotenpunktsverschiebungen des NodeContainers im Format:
        %
        %   [ u1 v1 w1;
        %     u2 v2 w2;
        %       ...
        %     un vn wn ]
        %
        function displacement_array = getDisplacementArray( self )
            number_of_nodes = self.getNumberOfNodes();
            dimension = self.node_list(1).dimension;
            
            % displacement_array is an ( number_of_nodes x dimension ) array.
            displacement_array = zeros( number_of_nodes, dimension );
            for i = 1:number_of_nodes
                displacement_array( i, : ) = self.node_list(i).dof.displacementFE;
            end
        end
        
        function setType( self, the_type )
            % This method not only sets the type but also initiates the list of all
            % soofea.model.RealIntegrationPoint objects. This is done by
            % iterating of the mathematical integration points soofea.numeric.NaturalIntegrationPoint inside
            % the class soofea.model.Type.
            
            % the_type: The actual soofea.model.type.Type of this object
            % (i.e. linear quad with two integration points per coordinate direction)
            self.type = the_type;
            natural_int_points = the_type.natural_int_points;
            
            for int_point = natural_int_points
                X = self.getCoordinateArray('material');                        % returns array with coordinates of node points
                H = self.type.shape.getArray( int_point.natural_coordinates );  % returns value of interpolation functions on the integration point on the reference cell
                coords = H*X;                                                   % calculates the coordinates of the RealIntegrationPoint with the isoparametric concept: x = H*X, just like u = H*û
                self.int_points = [self.int_points nsModel.RealIntegrationPoint( int_point, coords )];
            end
        end       
    end
end

