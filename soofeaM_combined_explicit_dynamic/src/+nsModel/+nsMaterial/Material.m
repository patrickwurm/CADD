classdef Material < nsModel.NumberedObject
    
    properties ( SetAccess = protected )
        two_dim_type
    end
    
    methods
        function self = Material(number, two_dim_type)
            self@nsModel.NumberedObject(number);

            if nargin < 2
                two_dim_type = {};
            end
            self.two_dim_type = two_dim_type;
 
        end
    end
end

