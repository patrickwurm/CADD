classdef TetraShape < nsModel.nsShape.Shape
    %Klasse für Triangle Shapes
        
    methods
        function self = TetraShape( order )
            % Die Dimension der Dreiecke ist 3!
            dimension = 3;
            
            self@nsModel.nsShape.Shape( order, dimension );
        end
    end
    
    methods ( Access = protected )
        % The value and derivative arrays are a bit more complicated for
        % tetras. See Hammer-Skriptum
        function val = getValue( self, natural_coordinates, node_index )
            I = node_index(1);
            J = node_index(2);
            K = node_index(3);
            L = ( self.order - ( node_index(1)-1 + node_index(2)-1 +node_index(3)-1 ) ) + 1;
            L_1 = natural_coordinates(1);
            L_2 = natural_coordinates(2);
            L_3 = natural_coordinates(3);
            L_4 = 1 - L_1 - L_2 - L_3;
            
            val = self.lagrange.getValue( L_1 , I , I ) *...
                  self.lagrange.getValue( L_2 , J , J ) *...
                  self.lagrange.getValue( L_3 , K , K ) *...
                  self.lagrange.getValue( L_4 , L , L );
        end
        
        % The value and derivative arrays are a bit more complicated for
        % tetras. See Hammer-Skriptum        
        function der = getDerivative( self, natural_coordinates, node_index, direction )
            L = ( self.order - ( node_index(1)-1 + node_index(2)-1 +node_index(3)-1 ) ) + 1;
            L_4 = 1-natural_coordinates(1)-natural_coordinates(2)-natural_coordinates(3);
            
            switch( direction )
                case 1
                    i = 1;
                    j = 2;
                    k = 3;
                case 2
                    i = 2;
                    j = 1;
                    k = 3;
                case 3
                    i = 3;
                    j = 1;
                    k = 2;
            end
            
            der = self.lagrange.getDerivative( natural_coordinates(i) , node_index(i) , node_index(i) ) * ...
                  self.lagrange.getValue( natural_coordinates(j) , node_index(j) , node_index(j) ) * ...
                  self.lagrange.getValue( natural_coordinates(k) , node_index(k) , node_index(k) ) * ...
                  self.lagrange.getValue( L_4 , L , L ) - ...
                  self.lagrange.getValue( natural_coordinates(i) , node_index(i) , node_index(i) ) * ...
                  self.lagrange.getValue( natural_coordinates(j) , node_index(j) , node_index(j) ) * ...
                  self.lagrange.getValue( natural_coordinates(k) , node_index(k) , node_index(k) ) * ...
                  self.lagrange.getDerivative( L_4 , L , L );
        end
        
        function node_positions = calcNodePositions(self)
            % Beim Tetraeder gehen die Koordinaten von 0 bis 1.
            nodes_per_side = self.order+1;
            node_positions = linspace(0,1,nodes_per_side);
        end
        
        function number_of_nodes = calcNumberOfNodes(self)
            % Die Gesamtanzahl der Knoten kann entsprechend folgender
            % Formel berechnet werden
            number_of_nodes = 2*(1 + self.order^2);
        end
        
        function nodeindex = getNodeIndex( self, local_node_number )
            if self.order == 1 
                %                 t  
                %               ,/
                %             ,/
                %            4                        
                %          ,/|`\                        
                %        ,/  |  `\              
                %      ,/    '.   `\                       
                %    ,/       |     `\             
                %  ,/         |       `\   
                % 1-----------'.--------3 --> s  
                %  `\.         |      ,/      
                %     `\.      |    ,/
                %        `\.   '. ,/
                %           `\. |/
                %              `2 
                %                 `\.
                %                    ` r
                switch local_node_number
                    case 1
                        nodeindex = [1,1,1];
                    case 2
                        nodeindex = [2,1,1];
                    case 3
                        nodeindex = [1,2,1];
                    case 4
                        nodeindex = [1,1,2];                        
                    otherwise
                        error('local node number out of range')
                end
            elseif self.order == 2
                %                 t  
                %               ,/
                %             ,/
                %            4                        
                %          ,/|`\                        
                %        ,/  |  `\              
                %      ,8    '.   `9                       
                %    ,/       |     `\             
                %  ,/         10       `\   
                % 1---------7-'.--------3 --> s  
                %  `\.         |      ,/      
                %     `\.      |    ,/
                %        5\.   '. ,6
                %           `\. |/
                %              `2 
                %                 `\.
                %                    ` r
                switch local_node_number
                    case 1
                        nodeindex = [1,1,1];
                    case 2
                        nodeindex = [3,1,1];
                    case 3
                        nodeindex = [1,3,1];
                    case 4
                        nodeindex = [1,1,3];
                    case 5
                        nodeindex = [2,1,1];
                    case 6
                        nodeindex = [2,2,1];
                    case 7
                        nodeindex = [1,2,1];
                    case 8
                        nodeindex = [1,1,2];
                    case 9
                        nodeindex = [1,2,2];
                    case 10
                        nodeindex = [2,1,2];                        
                    otherwise
                        error('local node number out of range')
                end
            else
                error('Tetras only implemented for order <= 2')
            end
        end
    end    
end