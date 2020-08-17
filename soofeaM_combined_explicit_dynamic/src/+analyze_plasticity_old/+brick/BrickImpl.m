classdef BrickImpl < analyze.ElementImpl
    methods
        function [k_element, f_int] = calcElementMatrices( self, element )
            global global_data;
            dim = global_data.dimension;
            
            n_amount = element.element_type.shape.node_amount;
            int_points = element.integration_points;
            
            %for i=1:n_amount
            %    int_points(i).disp();
            %end
            
            [int_tensor_1, int_tensor_2] = ...
                math.NumInt.integrateTensor( ...
                @(element,int_point)integrateConstitutiveComponent(self, element, int_point), ...
                element, int_points, [dim n_amount dim n_amount], [dim n_amount]);
            
            k_2d = tenmat(int_tensor_1, [1 2],[3 4] );
            k_element = k_2d.data;
            f_int = tenmat(int_tensor_2, 1, 2 );
        end
        
        function [int_tensor_1, int_tensor_2] = integrateConstitutiveComponent( self, element, int_point )
            global global_data;
            dim = global_data.dimension;
            
            %J = analyze.linear_tangent_map.LagrangianJacobian(element);
            %%lagrangian_j_inv has entries d r_k / d_xj
            %lagrangian_j_inv = J.getInv( int_point );
            %%der_tensor has entries d h_i/ d r_k
            der_tensor = element.element_type.shape.getDerivativeTensor(int_point.reference_coordinates);
            
            % Lagrangian jacobian
            node_point_tensor = element.getLagrangianCoord();
            J = ttt( node_point_tensor, der_tensor , 2 , 1 );
            jac_mat = tenmat( J, 1,2 );
            inverse_mat = inv(jac_mat.data);
            lagrangian_j_inv = tensor( inverse_mat );
            Det_J = det( jac_mat.data );
            if (Det_J<0)
                disp('Wrong local node numbering?');
                stopstop;
            end
            
            %ttt(A,B,2,1)
            %compute the contracted product of tensors A and B
            % Why 2,1 and not 2,2?
            % Second index of tensor A cancels with first index of tensor B!
            NJ = ttt(der_tensor,lagrangian_j_inv,2,1);
            
            delta = math.kronecker([dim dim]);
            delta_NJ = ttt(delta,NJ);
            
            buffer = 1/2 * ( permute(delta_NJ,[1 3 2 4]) + permute(delta_NJ,[1 3 4 2]) );
            
            delta_u = element.getEulerianCoord() - element.getLastConvergedCoord();
            delta_eps = ttt(buffer,delta_u,[1 2],[1 2]);
            
            u=element.getEulerianCoord() - element.getLagrangianCoord();
            eps = ttt(buffer,u,[1 2],[1 2]);
            
            props = [element.material.data_set.sigma_yield, element.material.data_set.iso_hardening_modulus ...
                 1/3*(3*element.material.data_set.bulk_modulus-2*element.material.data_set.shear_modulus), element.material.data_set.shear_modulus];
            stress_old = int_point.GetStressOld();
            state_old = int_point.GetStateOld();
            
            %CALL MATERIAL ROUTINE
            if element.material.data_set.material_type == 1
                [stress_new, state_new, dsde] = self.Umatel(delta_eps, props, stress_old, state_old);
            elseif element.material.data_set.material_type == 2
                [stress_new, state_new, dsde] = self.Umatiso(delta_eps, props, stress_old, state_old);
            else
                error(['Material #', num2str(element.material.data_set.material_type),' not implemented'])
            end
            
            int_point.SetStressNew(stress_new);
            int_point.SetStateNew(state_new);
            int_point.SetStrain(eps);
                
            c_buffer = ttt(dsde,buffer,[3 4],[3 4]);
            
            int_tensor_1 = ttt(buffer,c_buffer,[3 4],[1 2]) * Det_J;
            
            der_tensor_X = ttt( der_tensor, lagrangian_j_inv, 2, 1 );
            
            int_tensor_2 = ttt(stress_new, der_tensor_X, 2, 2) * Det_J;
            
            
        end
    end
    methods(Static)
        [stress_new, state_new, dsde] = Umatel(delta_eps, props, stress_old, state_old)
        [stress_new, state_new, dsde] = Umatiso(delta_eps, props, stress_old, state_old)
    end
end
