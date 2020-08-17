% ------------------------------------------------------------------
% This file is part of SOOFEAM -
%         Software for Object Oriented Finite Element Analysis in Matlab.
%
% SOOFEAM is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% SOOFEAM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with SOOFEAM.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------

classdef BCHandler < handle
    properties
        model;
        model_progression;
    end
    
    methods
        function self = BCHandler( model , model_progression )
            self.model = model;
            self.model_progression = model_progression;
        end
        
        function integrateBC( self , iteration_step, ...
                convergence_step, number_of_convergence_steps )
            %
            load_factor = convergence_step/number_of_convergence_steps;
            %
            prescribed_loads = self.model_progression.prescribed_load_cell;
            for i=1:length(prescribed_loads)
                self.integratePrescribedLoad( prescribed_loads{i}, load_factor );
            end
            %
            prescribed_dofs = self.model_progression.prescribed_dof_cell;
            if (iteration_step==1)
                load_factor = 1.0/number_of_convergence_steps;
            else
                load_factor = 0.0;
            end
            for i=1:length(prescribed_dofs)
                self.integratePrescribedDOF( prescribed_dofs{i}, load_factor );
            end
            %
        end
    end
    
    methods( Access = protected )
        function integratePrescribedDOF( self , prescribed_dof, load_factor )
            coord_ids = fieldnames(prescribed_dof.value_struct);
            for i=1:length(coord_ids)
                switch( coord_ids{i} )
                    case 'x'
                        self.model.getNode( prescribed_dof.node_number ). ...
                            getDisplacementDOF().setConstraintValue(1,prescribed_dof.value_struct.x * load_factor);
                    case 'y'
                        self.model.getNode( prescribed_dof.node_number ). ...
                            getDisplacementDOF().setConstraintValue(2,prescribed_dof.value_struct.y * load_factor);
                    case 'z'
                        self.model.getNode( prescribed_dof.node_number ). ...
                            getDisplacementDOF().setConstraintValue(3,prescribed_dof.value_struct.z * load_factor);                        
                end
            end
        end
        
        function integratePrescribedLoad( self , prescribed_load, load_factor )
            self.model.getNode(prescribed_load.node_number).setExtForce(prescribed_load.load_vector * load_factor);
        end
    end
end