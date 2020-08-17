classdef (Abstract) NumberedObject < handle
    
    properties ( SetAccess=private )
        number
    end
    
    methods
        function self = NumberedObject(number)
            self.number = number;
        end
    end
end

