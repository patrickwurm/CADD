classdef BoundaryImpl < handle
    properties
    end
    
    methods
        function surface_load = calcSurfaceLoad(self, boundary_component)
            surface_load = nsNumeric.NumInt.methodIntegrate( @self.surfaceLoadIntegrator, boundary_component.int_points, struct('boundary_component',boundary_component) );
        end
    end
    
    methods ( Static, Access = protected )
        function surface_force = surfaceLoadIntegrator( int_point, parameter )
            
            boundary_component = parameter.boundary_component;
            dimension = boundary_component.node_list(1).dimension;
            number_of_nodes = boundary_component.type.shape.getNumberOfNodes();
            n_dofs = number_of_nodes*dimension;
            load = int_point.surface_force;
            
            % calculate the force load only if the surface_force of the int_point is not empty.
            if ~isempty(load)
                % Get the relevant information for the boundary components integration points
                % H is a (1 x nodes_of_bc) array
                H = boundary_component.type.shape.getArray( int_point.getNaturalCoordinates() );
                                
                % jacobian is a (dimension_of_body x dimension_of_natural_boundary_element) array
                dV = nsAnalyzer.nsJacobian.BoundaryJacobian( boundary_component, int_point, 'material' ).getDet();                
                
                array = H' * load * dV;
                surface_force = reshape(array',[n_dofs,1]);
            else
                surface_force = zeros(n_dofs,1);
            end
        end
    end
end