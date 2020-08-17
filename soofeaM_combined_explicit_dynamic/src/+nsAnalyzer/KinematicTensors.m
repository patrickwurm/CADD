classdef KinematicTensors
    %KINEMATICTENSORS is a static-only class for the calculation of
    %important kinematic tensors. - i.e. more or less a namespace
    
    methods (Static)
        function F = calcDeformationGradient(material_jacobian, spatial_jacobian)
            j = transpose(spatial_jacobian.J_matrix);
            J = transpose(material_jacobian.J_matrix);
            F = j*inv(J); %#ok<MINV>
        end
        
        function C = calcRightCauchyGreen(deformation_gradient)
            F = deformation_gradient;
            C = transpose(F)*F;
        end
        
        function E = calcLagrangeStrain(right_cauchy_green)
            C = right_cauchy_green;
            dim = length(C);
            I = eye(dim);
            
            E = 0.5*( C - I );
        end
    end
end