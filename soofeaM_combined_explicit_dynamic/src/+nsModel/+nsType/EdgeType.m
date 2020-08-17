classdef EdgeType < nsModel.nsType.Type

    properties
    end
    
    methods
        function self = EdgeType(number, shape_order, number_of_int_points)
            self@nsModel.nsType.Type(number, shape_order, 'linear', number_of_int_points);
        end
    end
end

