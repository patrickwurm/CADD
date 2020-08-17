classdef ElementJacobian < nsAnalyzer.nsJacobian.Jacobian
    methods
        function self = ElementJacobian(node_container, int_point, configuration)
            self@nsAnalyzer.nsJacobian.Jacobian(node_container, int_point, configuration);
        end
        
        function Jinv = getInv(self)
            Jinv = inv(self.J_matrix);
        end
        
        function Det = getDet(self)
            Det = det(self.J_matrix);
        end
    end 
end