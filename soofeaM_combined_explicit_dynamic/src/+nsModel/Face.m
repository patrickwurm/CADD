classdef Face < nsModel.BoundaryComponent
    
    methods
        function self = Face( number, node_number_list )
            self@nsModel.BoundaryComponent( number, node_number_list );
        end
    end
end