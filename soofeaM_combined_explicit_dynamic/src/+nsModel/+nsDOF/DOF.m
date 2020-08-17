classdef DOF < handle
    % Degree Of Freedom Base Class - In multi-field formulations, there
    % might be more than just the #dimension displacement dofs. Derive
    % specific DOFs from this class.
    
    properties ( SetAccess = protected )
        constraint
    end
    
    methods ( Abstract )
        getConstraint(self)
        resetDOF(self)
    end    
end

