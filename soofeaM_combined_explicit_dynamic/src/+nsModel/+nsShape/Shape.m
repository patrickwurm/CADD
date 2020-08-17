classdef Shape < handle
    % Abstrakte Basisklasse für Ansatzfunktionen
    
    properties ( SetAccess = protected )
        order
        dimension
        lagrange
    end
    
    methods
        function self = Shape(order, dimension)
            self.order = order;
            self.dimension = dimension;
            self.lagrange = nsModel.nsShape.Lagrange(self.calcNodePositions());
        end
              
        % The H-matrix or interpolation matrix has the form
        %
        % H = [h1 h2 ... hn]
        %
        % where n is the number of nodes and hi is the i-th Ansatzfunction
        function H = getArray(self, coordinates)
            N = self.calcNumberOfNodes;
            H = zeros(1,N);
            for node_number=1:N
                node_index = self.getNodeIndex(node_number);
                H(1,node_number) = self.getValue( coordinates, node_index  );
            end
        end
              
        % The dH-matrix of derivative matrix has the form
        %
        % dH = [ dh1_dr dh2_dr ... dhn_dr
        %        dh1_ds dh2_ds ... dhn_ds
        %       (dh1_dt dh2_dt ... dhn_dt) ]
        %
        % where n is the number of nodes and dhi_drj  is the derivative of
        % the i-th Ansatzfunction in the rj coordinate direction.
        %
        % IMPORTANT: The derivatives are with respect to the natural
        % coordinates! In order to calculate stiffness matrices, we will
        % need jacobians.
        function dH = getDerivativeArray(self, natural_coordinates)
            N = self.calcNumberOfNodes;
            dH = zeros(self.dimension, N);
            for node_number=1:N
                node_index = self.getNodeIndex(node_number);
                for direction=1:self.dimension
                    dH(direction, node_number) = self.getDerivative( natural_coordinates, node_index, direction );
                end
            end
        end
        
        function number_of_nodes = getNumberOfNodes(self)
           number_of_nodes = self.calcNumberOfNodes;
        end
    end    
    
    methods (Abstract, Access = protected)
        getValue(self, natural_coordinates, node_index)
        getDerivative(self, natural_coordinates, node_index, derivative_direction)
        calcNodePositions(self)
        calcNumberOfNodes(self)
        getNodeIndex(local_node_number)
    end
end