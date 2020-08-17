classdef FaceType < nsModel.nsType.Type
    
    properties
    end
    
    methods
        function self = FaceType(number, shape_order, shape_type, number_of_int_points)
            self@nsModel.nsType.Type(number, shape_order, shape_type, number_of_int_points);
        end
    end
end

