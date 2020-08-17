classdef ElementImpl
    methods
        function F_body = calcBodyForce(self, element)
            % Body forces are integrated over the volume of the body.
            F_body = nsNumeric.NumInt.methodIntegrate(@self.bodyForcesIntegrator, element.int_points, struct('element',element));
        end
    
        % for stress projection
        function M = calcMassMatrix(self, element, params)
            if nargin < 3 %used for stress projection (only 1 variable is interpolated)
                M = nsNumeric.NumInt.methodIntegrate( @self.massMatrixIntegrator, element.int_points, struct('element',element) );
            else
                M = nsNumeric.NumInt.methodIntegrate( @self.massMatrixIntegrator, element.int_points, struct('element',element,'full_rank',params) );
            end
        end
        
        % for stress projection
        function proj_vec = calcProjectionVector(self, element)
            proj_vec = nsNumeric.NumInt.methodIntegrate( @self.projectionVectorIntegrator, element.int_points, struct('element',element) );
        end
    end
       
    methods (Static)
        % calc comparison stresses in IP
        function calcComparisonStressesInIP(element)
            for int_point = element.int_points
                
                stress_3D = zeros(3,3);
                
                if element.node_list(1).dimension == 2
                    if strcmp(element.material.two_dim_type,'plane_stress')
                        stress_3D(1,1) = int_point.stress(1);
                        stress_3D(2,2) = int_point.stress(2);
                        stress_3D(1,2) = int_point.stress(3);
                        stress_3D(2,1) = stress_3D(1,2);%other entries remain zero
                    elseif strcmp(element.material.two_dim_type,'plane_strain')
                        stress_3D(1,1) = int_point.stress(1);
                        stress_3D(2,2) = int_point.stress(2);
                        stress_3D(1,2) = int_point.stress(3);
                        stress_3D(2,1) = stress_3D(1,2);
                        stress_3D(3,3) = element.material.calcThirdNormalStressComponentForPlaneStrain( int_point );
                    end
                elseif element.node_list(1).dimension == 3
                    stress_3D(1,1) = int_point.stress(1);
                    stress_3D(2,2) = int_point.stress(2);
                    stress_3D(3,3) = int_point.stress(3);
                    stress_3D(1,2) = int_point.stress(4);
                    stress_3D(2,3) = int_point.stress(5);
                    stress_3D(1,3) = int_point.stress(6);
                    stress_3D(2,1) = stress_3D(1,2);
                    stress_3D(3,2) = stress_3D(2,3);
                    stress_3D(3,1) = stress_3D(1,3);
                else
                    error('Cauchy stress evaluation only implemented for 2D and 3D problems')
                end
                
                principal_stresses = eig(stress_3D);
                
                int_point.setStressVonMises( sqrt(1/2*( (principal_stresses(1)-principal_stresses(2))^2 +...
                    (principal_stresses(2)-principal_stresses(3))^2 +...
                    (principal_stresses(3)-principal_stresses(1))^2 )) );
            end
        end
    end
    
    methods (Access = protected, Static)      
        function M = massMatrixIntegrator( int_point, parameters )
            element = parameters.element;
            H = element.type.shape.getArray( int_point.getNaturalCoordinates() );
            
            if isfield(parameters,'full_rank')
                blocks = cell(1,length(H));
                for i=1:length(H) 
                blocks{i} = H(i)*eye(element.type.shape.dimension);
                end
                H = [];
                for i=1:numel(blocks)
                    H = horzcat(H, blocks{i});
                end
            end
            M = H' * H * int_point.dV;
        end
        
        function proj_vec = projectionVectorIntegrator(int_point, parameters )
            element = parameters.element;
            H = element.type.shape.getArray( int_point.getNaturalCoordinates() );

            proj_vec = H' * int_point.stress_von_mises * int_point.dV;
        end
        
        function body_force = bodyForcesIntegrator( int_point, parameter )
            
            element = parameter.element;
            dimension = element.node_list(1).dimension;
            number_of_nodes = element.type.shape.getNumberOfNodes();
            n_dofs = number_of_nodes*dimension;
            load = int_point.body_force;
            
            % calculate the force load only if the surface_force of the int_point is not empty.            
            if ~isempty(load)
                % Get the relevant information for the elements integration points
                % H is a (1 x nodes_of_element) array
                H = element.type.shape.getArray( int_point.getNaturalCoordinates() );
                
                % jacobian is a (dimension_of_body x dimension_of_element)
                % array (hidden in int_point.dV)
                array = H' * load *  int_point.dV;
                body_force = reshape(array',[n_dofs, 1]);
            else
                body_force = zeros(number_of_nodes*dimension,1);
            end
        end
    end
end
