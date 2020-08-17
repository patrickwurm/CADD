classdef PadAtom < nsModel.NumberedObject
    
    properties( SetAccess = protected )
        coordinates
        natural_coordinates
        displacements
    end
    
    methods
        function self = PadAtom( number, coordinates )
            self@nsModel.NumberedObject(number);
            self.coordinates = coordinates;
            self.natural_coordinates = [];
        end
             
        function setNaturalCoords( self, coords)
            self.natural_coordinates = coords;
        end
        
        function setDisplacements( self, displacements)
            self.displacements = displacements;
        end
        
    end
end