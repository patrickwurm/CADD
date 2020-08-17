classdef Type < nsModel.NumberedObject
    
    properties ( SetAccess = protected )
        % The implementation gets set in the concrete example file.
        implementation
        
        % Shape and natural integration points get set in the constructor
        % below.
        shape
        natural_int_points
    end
    
    methods
        % The type-object is created in the concrete example file.
        function self = Type(number, shape_order, shape_type, number_of_int_points)
            self@nsModel.NumberedObject(number);
            
            % Set the shape.
            if strcmp(shape_type, 'linear')
                self.shape = nsModel.nsShape.LinearShape(shape_order);
            elseif strcmp(shape_type,'quad')
                self.shape = nsModel.nsShape.QuadShape(shape_order);
            elseif strcmp(shape_type, 'triangle')
                self.shape = nsModel.nsShape.TriangleShape(shape_order);
            elseif strcmp(shape_type, 'hex')
                self.shape = nsModel.nsShape.HexShape(shape_order);
            elseif strcmp(shape_type, 'tetra')
                self.shape = nsModel.nsShape.TetraShape(shape_order);
            else
                error(['The given shape type ',shape_type,' is not implemented'])
            end
            
            % Create the natural integration points on the reference cell
            % for this type.
            self.natural_int_points = nsNumeric.NumInt.getNaturalIntegrationPoints(shape_type, number_of_int_points);
        end       
        
        function setImplementation(self, implementation)
            self.implementation = implementation;
        end
    end    
end