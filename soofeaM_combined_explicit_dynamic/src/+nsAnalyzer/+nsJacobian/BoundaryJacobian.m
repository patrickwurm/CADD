classdef BoundaryJacobian < nsAnalyzer.nsJacobian.Jacobian
    properties (SetAccess = protected)
        tangent1
        tangent2
    end
    
    methods
        function self = BoundaryJacobian(node_container, int_point, configuration)
            self@nsAnalyzer.nsJacobian.Jacobian(node_container, int_point, configuration);
            self.calcTangentVectors();
        end
             
        function Det = getDet(self)
            % Die Determinante wird benötigt, um das Integral auf das
            % Einheitselement zu transformieren.
            J = self.J_matrix;
            dimJ = size(J);
            if min(dimJ) == 2
                % In diesem Fall beschreibt die Jacobi-Matrix die
                % Transformation eines Faces: (x,y,z) -> (r,s)
                Det = norm(cross(J(1,:), J(2,:)));
            elseif min(dimJ) == 1
                % In diesem Fall beschreibt die Jacobi-Matrix die
                % Transformation eines Edges: (x,y) -> (r)
                Det = norm(J(1,:));
            else
                error('Error in Jacobian-dimensions!');
            end
        end
        
        function normal = getNormal(self)
            normal = cross(self.tangent1, self.tangent2);
        end
    end
    
    methods (Access = private)
        function calcTangentVectors(self)
            dimJ = size(self.J_matrix);
            
             if min(dimJ) == 2
                % In diesem Fall beschreibt die Jacobi-Matrix die
                % Transformation eines Faces: (x,y,z) -> (r,s)
                self.tangent1 = self.J_matrix(1,:);
                self.tangent2 = self.J_matrix(2,:);
            elseif min(dimJ) == 1
                % In diesem Fall beschreibt die Jacobi-Matrix die
                % Transformation eines Edges: (x,y) -> (r)
                self.tangent1 = self.J_matrix(1,:);
                self.tangent2 = [0;0;1];
             end
        end
    end 
end