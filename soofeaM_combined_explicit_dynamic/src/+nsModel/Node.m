classdef Node < nsModel.NumberedObject
    %NODE Node of the finite element grid
    
    properties
        % needed for a nonlinear formulation.
        % still there for the linear one.
        material_coordinates
        spatial_coordinates
        
        % needed for a correct plastic implementation.
        last_converged_coordinates
        
        % formulation-independent properties
        dimension
        dof
        stress_von_mises
        node_force
        interface_atom_number
    end
    
    methods
        % A node gets created only with its material coordinates, at the
        % initial state t=0.
        function self = Node( number, coordinates )
            self@nsModel.NumberedObject(number);
            self.material_coordinates = coordinates;
            self.spatial_coordinates = coordinates;
            self.last_converged_coordinates = coordinates;
            self.dimension = length(coordinates);
            
            self.dof = nsModel.nsDOF.DisplacementDOF(self.dimension);
            
            % Wird vorerst für alle Knoten zu 0 gesetzt. Später werden
            % eventuelle Node Forces gesetzt.
            self.node_force  = zeros(1,self.dimension);
            
            % Zu Beginn wird von einem spannungsfreien Zustand ausgegangen
            self.stress_von_mises = 0;
        end
        
        % call this method in the BCHandler of a linear analysis where
        % absolute displacements are Dirichlet boundary conditions
        function setBCDisplacement( self, coord_flag, displacement )
            % Pass a coordinate flag ('x', 'y' or 'z') to decide in which
            % coordinate direction the displacement is constrained.
            switch coord_flag
                case 'x'
                    self.dof.setConstraintDisplacement( 1, displacement );
                case 'y'
                    assert(self.dimension>1,'Prescribe y-value only if dimension of problem is > 1!');
                    self.dof.setConstraintDisplacement( 2, displacement );
                case 'z'
                    assert(self.dimension>2,'Prescribe z-value only if dimension of problem is > 2!');
                    self.dof.setConstraintDisplacement( 3, displacement );
            end
        end
        
        % call this method in the BCHandler of a NONlinear analysis where
        % incremental displacements are Dirichlet boundary conditions
        function setBCIncrement( self, coord_flag, increment )
            % Pass a coordinate flag ('x', 'y' or 'z') to decide in which
            % coordinate direction the increment is constrained.
            switch coord_flag
                case 'x'
                    self.dof.setConstraintIncrement( 1, increment );
                case 'y'
                    assert(self.dimension>1,'Prescribe y-value only if dimension of problem is > 1!');
                    self.dof.setConstraintIncrement( 2, increment );
                case 'z'
                    assert(self.dimension>2,'Prescribe z-value only if dimension of problem is > 2!');
                    self.dof.setConstraintIncrement( 3, increment );
            end
        end
           
        % needed for the nonlinear analyis. called once per newton
        % iteration. call this method ONLY in the nonlinear analysis!
        % This method is called once per newton-step
        function updateCoordinates( self )
            self.spatial_coordinates = self.material_coordinates + self.dof.displacement;
            % alternative:
            % self.spatial_coordinates += self.dof.increment
        end        
        
        % Node Forces are set in the BCHandler
        function setNodeForce( self, node_force )
            self.node_force = node_force;
        end
        
        function setStressVonMises(self, svm)
            self.stress_von_mises = svm;
        end
        
        function setInterfaceAtomNumber( self, atom_number )
            self.interface_atom_number = atom_number;
        end
    end
end

