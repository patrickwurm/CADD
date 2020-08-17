classdef Jacobian < handle
    properties ( SetAccess = protected )
        node_container
        node_coordinates
        int_point
        J_matrix
    end
    
    methods
        function self = Jacobian(node_container, int_point, configuration)
            self.node_container = node_container;
            self.int_point = int_point;
            self.node_coordinates = node_container.getCoordinateArray(configuration);
            self.calc()
        end
    end
    
    methods ( Access = private )
        function calc(self)
        % Diese Methode wird vom Konstruktor aus aufgerufen und berechnet
        % die Jacobimatrix im Integrationspunkt.
            dH = self.node_container.type.shape.getDerivativeArray(self.int_point.getNaturalCoordinates());
            self.J_matrix = dH*self.node_coordinates;
        end
    end
    
   methods ( Static )
       function [J_matrix] = staticCalcJacobian(node_container, int_point, configuration)
           node_coordinates = node_container.getCoordinateArray(configuration);
           J_matrix = int_point.der_array*node_coordinates;
       end
   end
end