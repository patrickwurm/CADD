classdef DisplacementDOF < nsModel.nsDOF.DOF
    % specifies displacement DOFs
    
    % NAME IS MISLEADING (maybe?) CAUSE THERE ARE ALSO VELOCITIES AND
    % ACCELERATIONS INVOLVED
    
    properties ( SetAccess = protected )
        displacement
        displacementFE
        displacementincrement
        velocity
        acceleration
    end
    
    methods
        function self = DisplacementDOF( dimension )
            % the values of the displacementDOFs are the displacements. they are stored in an 1xdim array.
            self.displacement = zeros( 1, dimension );
            % the increments are 0, too, at the beginning.
            % note: the increment is the tiny increment of the current
            % newton-iteration, not all cumulated increments of the current
            % time step!
            self.displacementincrement = zeros( 1, dimension );
            % Per default, all values of a displacement DOF are not constraint,
            % hence they are initiated as 'false'.
            self.displacementFE = zeros( 1, dimension );
            self.velocity = zeros( 1, dimension );
            self.acceleration = zeros( 1, dimension );
            self.constraint = false( 1, dimension );
        end
        
        function displacement = getDisplacement( self, coordinate_id )
            displacement = self.displacement( coordinate_id );
        end
        
        function inc_disp = getIncrementalDisplacement(self, coordinate_id)
            inc_disp = self.displacementincrement( coordinate_id );
        end
        
        function setDisplacement(self, coord_id, displacement)
           self.displacement(coord_id) = displacement; 
        end 
        
        function constraint = getConstraint( self, coordinate_id )
            constraint = self.constraint( coordinate_id );
        end
              
        function setConstraintDisplacement( self, coordinate_id, displacement )
            self.displacement (coordinate_id ) = displacement;
            self.constraint( coordinate_id ) = true;
        end
        
        function setConstraintIncrement( self, coordinate_id, increment )
            self.displacementincrement (coordinate_id ) = increment;
            self.constraint( coordinate_id ) = true;
        end
        
        function resetDOF( self, coordinate_id )
            self.displacementincrement(coordinate_id) = 0;
            self.constraint(coordinate_id) = false;
        end
        
        function addIncrement( self, coordinate_id, increment )
            self.displacement(coordinate_id) = self.displacement(coordinate_id) + increment;
            self.displacementincrement(coordinate_id) = increment;
        end

        function velocity = getVelocity( self, coordinate_id )
            velocity = self.velocity( coordinate_id );
        end
        
        function acceleration = getAcceleration( self, coordinate_id )
            acceleration = self.acceleration( coordinate_id );
        end
        
        function setVelocity(self, coord_id, velocity)
           self.velocity(coord_id) = velocity; 
        end
        
        function setAcceleration(self, coord_id, acceleration)
            self.acceleration(coord_id) = acceleration;
        end
        
        function setDisplacementFE(self, coord_id, displacement)
            self.displacementFE(coord_id) = displacement;
        end
        
    end
end