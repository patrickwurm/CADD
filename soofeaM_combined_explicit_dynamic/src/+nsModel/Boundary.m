classdef Boundary < nsModel.NodeContainer
    
    properties ( SetAccess = private )
        % This list contains all boundary components. This might be
        % soofea.model.Point (currently not implemented) for the 1D case,
        % soofea.model.Edge for the 2D case or soofea.model.Face for the 3D
        % case.
        component_list
    end
    
    methods
        function self = Boundary( number )
            % The following lines are neccessary, because the boundaries
            % may be constructed unordered, e.g. 1-3-2-..
            % see: https://de.mathworks.com/help/matlab/matlab_oop/initialize-object-arrays.html
            if nargin == 0
                number = 0;
            end
            % The node_number_list is empty at the point of creation of the
            % boundary
            self@nsModel.NodeContainer( number, [] );
        end
        
        function addComponent( self, new_component )
            % As this class is derived from soofea.model.NodeContainer, we have
            % to append the nodes of the new component and their numbers to the
            % nsModel.NodeContainer.node_number_list and the
            % nsModel.NodeContainer.node_list.
            
            % The component given to this function (soofea.model.Edge in 2D,
            % soofea.model.Face in 3D) is added to the boundary
            self.component_list = [self.component_list new_component];
            
            for new_node = new_component.node_list              
                % Add the new node only if it is not already contained and
                % if the current new_node is not empty.
                % (Why should it be empty?)
                if ~ismember(new_node, self.node_list) && ~isempty(new_node) 
                    self.node_list        = [self.node_list new_node];
                    self.node_number_list = [self.node_number_list new_node.number];
                end
            end
        end
    end
end