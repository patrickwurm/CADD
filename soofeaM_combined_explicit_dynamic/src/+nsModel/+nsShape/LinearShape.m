classdef LinearShape < nsModel.nsShape.Shape
    % Klasse für eindimensionale Ansatzfunktionen im Referenzintervall
    
    properties
    end
    
    methods
        function self = LinearShape(order)
            % Bei 1D-Ansatzfunktionen ist die Dimension natürlich 1.
            dimension = 1;
           
            self@nsModel.nsShape.Shape(order, dimension)
        end
    end
    
    methods ( Access = protected )
        function val = getValue(self, natural_coordinate, node_index)
            % Der Funktionswert entspricht dem Wert des Lagrange-Polynoms
            val = self.lagrange.getValue(natural_coordinate, node_index);
        end
        
        function der = getDerivative(self, natural_coordinate, node_index, ~)
            % Die Ableitung entspricht der des Lagrange-Polynoms. Der
            % letzte Funktionsparameter - die Richtung der Ableitung - wird
            % im 1D-Fall nicht benötigt, da es nur eine Richtung gibt.
            der = self.lagrange.getDerivative(natural_coordinate, node_index);
        end
        
        function node_positions = calcNodePositions(self)
            % Die Stützstellen sind äquidistant über das Intervall [-1 1]
            % verteilt.
            n = self.calcNumberOfNodes;
            node_positions = linspace(-1,1,n);
        end
        
        function number_of_nodes = calcNumberOfNodes(self)
            % Im 1D-Fall ist die Anzahl der Stützstellen gleich wie bei
            % Lagrange-Polynomen:
            number_of_nodes = self.order + 1;
        end

        function nodeindex = getNodeIndex(self, local_node_number)
            % Der Index muss für Elemente höherer Ordnung ebenfalls
            % entsprechend der außen-innen-Regel berechnet werden.
            if self.order == 1
                nodeindex = local_node_number;
            elseif self.order == 2
                switch local_node_number
                    case 1
                        nodeindex = 1;
                    case 2
                        nodeindex = 3;
                    case 3
                        nodeindex = 2;
                end
            else
                error('LinearShape only implemented for order <= 2')
            end
        end
    end
end
