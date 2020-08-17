classdef BCHandler < handle
    
    properties (SetAccess = protected)
        model
    end
    
    methods
        function setModel(self, model)
            self.model = model;
        end
        
        % at the nonlinear analysis: resets at the beginning of each time
        % step: increment is set to '0', constraint is set to 'false'.
        function resetBC(self)
            for node = self.model.node_dict
                for coord_id = 1:self.model.dimension
                    node.dof.resetDOF(coord_id);
                end
            end
        end
        
        % done after the first newton-iteration -> once per time-step
        % after the very first newton step, the boundary
        % conditions are set to zero! explain better...
        function setPrescribedDOFZero(self)
            for node = self.model.node_dict
                for coord_id = 1:self.model.dimension
                    if node.dof.getConstraint(coord_id)
                        node.dof.setConstraintIncrement(coord_id, 0.0)
                    end
                end
            end
        end
    end
    
    % This method is implemented for each example seperately, and
    % specifies the boundary conditions for the specific geometry and load
    % case.
    methods (Static, Abstract)
        incorporateBC(model, time, old_matrix)
    end

end