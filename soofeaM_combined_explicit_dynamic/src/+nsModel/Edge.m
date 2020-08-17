classdef Edge < nsModel.BoundaryComponent
    
    methods
        function self = Edge( number, node_number_list )
            self@nsModel.BoundaryComponent( number, node_number_list );
        end
    end
end