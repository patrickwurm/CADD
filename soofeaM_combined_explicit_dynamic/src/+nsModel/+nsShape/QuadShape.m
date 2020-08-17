classdef QuadShape < nsModel.nsShape.Shape
    % Klasse für 2d-Ansatzfunktionen am Referenzviereck.
    
    methods
        function self = QuadShape(order)
            % Die Dimension der Rechteckelemente ist 2.
            dimension = 2;
            
            self@nsModel.nsShape.Shape(order, dimension)
        end
    end
    
    methods ( Access = protected )
        function val = getValue(self, natural_coordinates, node_index)
            % Der Funktionswert entspricht der Theorie:
            % h[i,j](r,s) = h_i(r)*h_j(s)
            val = self.lagrange.getValue(natural_coordinates(1), node_index(1)) * ...
                self.lagrange.getValue(natural_coordinates(2), node_index(2));
        end
        
        function der = getDerivative(self, natural_coordinates, node_index, derivative_direction)
            % Die Ableitung entspricht ebenfalls der Theorie.
            % Zusätzlich muss die Richtung der Ableitung angegeben werden:
            % Ableitung nach r->1 oder nach s->2
            if derivative_direction == 1
                der = self.lagrange.getDerivative(natural_coordinates(1), node_index(1)) * ...
                    self.lagrange.getValue(natural_coordinates(2), node_index(2));
            elseif derivative_direction == 2
                der = self.lagrange.getValue(natural_coordinates(1), node_index(1)) * ...
                    self.lagrange.getDerivative(natural_coordinates(2), node_index(2));
            end
        end
        
        function node_positions = calcNodePositions(self)
            % Nachdem in beide Koordinatenrichtungen gleich viele
            % Knotenpunkte existieren, genügt es die Positionen in
            % eindimensionaler Richtung zu berechnen.
            nodes_per_side = self.order+1;
            node_positions = linspace(-1,1,nodes_per_side);
        end
        
        function number_of_nodes = calcNumberOfNodes(self)
            % Die Gesamtanzahl der Knoten kann entsprechend folgender
            % Formel berechnet werden
            number_of_nodes = (self.order + 1)^2;
        end
        
        function node_index = getNodeIndex(self, local_node_number)
            if self.order == 1
                %        s
                %       ^
                %       |
                % 4-----------3
                % |     |     |
                % |     |     |
                % |     +---- | --> r
                % |           |
                % |           |
                % 1-----------2
                switch local_node_number
                    case 1
                        node_index = [1,1];
                    case 2
                        node_index = [2,1];
                    case 3
                        node_index = [2,2];
                    case 4
                        node_index = [1,2];
                    otherwise
                        error('local node number out of range')
                end
            elseif self.order == 2
                %        s
                %       ^
                %       |
                % 4-----7-----3
                % |     |     |
                % |     |     |
                % 8     9---- 6 --> r
                % |           |
                % |           |
                % 1-----5-----2
                switch local_node_number
                    case 1
                        node_index = [1,1];
                    case 2
                        node_index = [3,1];
                    case 3
                        node_index = [3,3];
                    case 4
                        node_index = [1,3];
                    case 5
                        node_index = [2,1];
                    case 6
                        node_index = [3,2];
                    case 7
                        node_index = [2,3];
                    case 8
                        node_index = [1,2];
                    case 9
                        node_index = [2,2];
                    otherwise
                        error('local node number out of range')
                end
            else
                error('QuadShape only implemented for order <= 2')
            end
        end
    end
end
