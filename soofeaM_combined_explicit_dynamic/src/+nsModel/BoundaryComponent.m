classdef BoundaryComponent < nsModel.NodeContainer    
    methods
        function self = BoundaryComponent( number, node_number_list )
            self@nsModel.NodeContainer( number, node_number_list )
        end
    end    
end

