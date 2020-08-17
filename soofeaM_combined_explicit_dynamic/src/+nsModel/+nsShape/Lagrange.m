classdef Lagrange < handle
    %LAGRANGE Klasse zur Implementierung von eindimensionalen Lagrange-Polynomen.
    
    properties ( SetAccess = protected )
        node_positions
        order
    end
    
    methods
        function self=Lagrange(node_positions)
            % node_positions: Liste der Stützstellen des Polynoms
            self.node_positions = node_positions;
            self.order = length(self.node_positions) - 1;
        end
        
        function L = getValue(self, coordinate, index, top_index)
            % Gibt den Wert des Lagrangepolynoms mit Index 'index' an der Stelle 'r=coordinate' zurück.
            % Der Wert top_index wird für allgemeine Dreieckselemente
            % benötigt.
            
            if nargin <= 3 % nargin == 2 wird für lineare, quad und hex Elemente benötigt.
                high = self.order + 1;
            else %nargin == 3 wird für Triangle- und Tetraelementen benötigt. Bei diesen Elementen
                 %variiert der Grad des Lagrangepolynoms je nach Richtung und Knoten und ist somit nicht
                 %immer gleich self.order+1;
                high = top_index;
            end
            
            L = 1;            
            for i=1:high
                if i~=index
                    L = L*(coordinate-self.node_positions(i))/(self.node_positions(index)-self.node_positions(i));
                end
            end
        end
        
        function D = getDerivative(self, coordinate, index, top_index)
            % Gibt den Wert der Ableitung des Lagrangepolynoms mit Index 'index' an der Stelle 'r=coordinate' zurück.
            % Der Wert top_index wird für allgemeine Dreieckselemente
            % benötigt.
            
            if nargin <= 3 % nargin == 2 wird für lineare, quad und hex Elemente benötigt.
                high = self.order + 1;
            else %nargin == 3 wird für Triangle- und Tetraelementen benötigt. Bei diesen Elementen
                 %variiert der Grad des Lagrangepolynoms je nach Richtung und Knoten und ist somit nicht
                 %immer gleich self.order+1;
                high = top_index;
            end

            factor_1 = 1;
            for i=1:high
                if i~=index
                    factor_1 = factor_1/(self.node_positions(index) - self.node_positions(i));
                end
            end
            sum_i = 0;
            
            for i=1:high
                if i~=index
                    inner_factor = 1;
                    for j=1:self.order+1
                        if j~=index && j~=i
                            inner_factor = inner_factor * (coordinate - self.node_positions(j));
                        end
                    end
                    sum_i = sum_i + inner_factor;
                end
            end
            D = factor_1 * sum_i;
        end
    end
end