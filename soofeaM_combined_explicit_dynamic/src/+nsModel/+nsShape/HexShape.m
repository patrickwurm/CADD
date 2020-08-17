classdef HexShape < nsModel.nsShape.Shape
    % Klasse für 3d-Ansatzfunktionen am Referenzwürfel.
    
    methods
        function self = HexShape(order)
            % Die Dimension der Würfelelemente ist 3.
            dimension = 3;
            
            self@nsModel.nsShape.Shape(order, dimension)
        end
    end
    
    methods ( Access = protected )
        function val = getValue(self, natural_coordinates, node_index)
            % Der Funktionswert entspricht der Theorie:
            % h[i,j,k](r,s,t) = h_i(r)*h_j(s)*h_k(t)
            val = self.lagrange.getValue(natural_coordinates(1), node_index(1)) * ...
                self.lagrange.getValue(natural_coordinates(2), node_index(2)) * ...
                self.lagrange.getValue(natural_coordinates(3), node_index(3));
        end
        
        function der = getDerivative(self, natural_coordinates, node_index, derivative_direction)
            % Die Ableitung entspricht ebenfalls der Theorie.
            % Zusätzlich muss die Richtung der Ableitung angegeben werden:
            % Ableitung nach r->1, nach s->2 oder nach t->3
            if derivative_direction == 1
                der = self.lagrange.getDerivative(natural_coordinates(1), node_index(1)) * ...
                    self.lagrange.getValue(natural_coordinates(2), node_index(2)) * ...
                    self.lagrange.getValue(natural_coordinates(3), node_index(3));
            elseif derivative_direction == 2
                der = self.lagrange.getValue(natural_coordinates(1), node_index(1)) * ...
                    self.lagrange.getDerivative(natural_coordinates(2), node_index(2)) * ...
                    self.lagrange.getValue(natural_coordinates(3), node_index(3));
            elseif derivative_direction == 3
                der = self.lagrange.getValue(natural_coordinates(1), node_index(1)) * ...
                    self.lagrange.getValue(natural_coordinates(2), node_index(2)) * ...
                    self.lagrange.getDerivative(natural_coordinates(3), node_index(3));
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
            number_of_nodes = (self.order + 1)^3;
        end
        
        function node_index = getNodeIndex(self, local_node_number)
            if self.order == 1
                %        t
                % 8----------7
                % |\  s  ^   |\
                % | \ \  |   | \
                % |  \ \ |   |  \
                % |   5-\----+---6
                % |   |  +-- |-- | -> r
                % 4---+------3   |
                %  \  |       \  |
                %   \ |        \ |
                %    \|         \|
                %     1----------2
                
                switch( local_node_number )
                    case 1
                        node_index = [1,1,1];
                    case 2
                        node_index = [2,1,1];
                    case 3
                        node_index = [2,2,1];
                    case 4
                        node_index = [1,2,1];
                    case 5
                        node_index = [1,1,2];
                    case 6
                        node_index = [2,1,2];
                    case 7
                        node_index = [2,2,2];
                    case 8
                        node_index = [1,2,2];
                    otherwise
                        error('local node number out of range')
                end
            elseif self.order == 2
            % needs to be tested! 
                % 8----20----7     
                % |\         |\         t
                % |18    26  | 19    s  |
                % 16 \ 25    15 \     \ |
                % |   5----17+---6     \|
                % |23 |  27  | 24|      +-----r
                % 4---+-14---3   | 
                %  \ 11    22 \  13
                %  10 |  21   12 | 
                %    \|         \| 
                %     1-----9----2 
                switch local_node_number
                    case 1
                        node_index = [1,1,1];
                    case 2
                        node_index = [3,1,1];
                    case 3
                        node_index = [3,3,1];
                    case 4
                        node_index = [1,3,1];
                    case 5
                        node_index = [1,1,3];
                    case 6
                        node_index = [3,1,3];
                    case 7
                        node_index = [3,3,3];
                    case 8
                        node_index = [1,3,3];
                    case 9
                        node_index = [2,1,1];
                    case 10
                        node_index = [1,2,1];
                    case 11
                        node_index = [1,1,2];
                    case 12
                        node_index = [3,2,1];
                    case 13
                        node_index = [3,1,2];
                    case 14
                        node_index = [2,3,1];
                    case 15
                        node_index = [3,3,2];
                    case 16
                        node_index = [1,3,2];
                    case 17
                        node_index = [2,1,3];
                    case 18
                        node_index = [1,2,3];
                    case 19
                        node_index = [3,2,3];
                    case 20
                        node_index = [2,3,3];
                    case 21
                        node_index = [2,2,1];
                    case 22
                        node_index = [2,1,2];
                    case 23
                        node_index = [1,2,2];
                    case 24
                        node_index = [3,2,2];
                    case 25
                        node_index = [2,3,2];
                    case 26
                        node_index = [2,2,3];
                    case 27
                        node_index = [2,2,2];
                    otherwise
                        error('local node number out of range')
                end
            else
                error('QuadShape only implemented for order <= 2')
            end
        end
    end
end
