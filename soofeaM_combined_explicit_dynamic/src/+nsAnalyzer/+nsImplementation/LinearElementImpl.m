classdef LinearElementImpl < nsAnalyzer.nsImplementation.ElementImpl
    properties
    end
    
    methods
        function calcStrainStressInIP(self, element)
            for int_point = element.int_points
                %jacobian = nsAnalyzer.nsJacobian.ElementJacobian(element, int_point, 'material');
                %jacobian = nsAnalyzer.nsJacobian.Jacobian.staticCalcJacobian(element, int_point, 'material')
                
                %int_point.dV = det(int_point.jacobian);
                %int_point.B = self.getBMatrix(int_point.jacobian, int_point.der_array, element);
                
                strain_voigt = int_point.B * self.coordinatesToVoigt(element.getDisplacementArray());
                int_point.setStress(int_point.CC * strain_voigt);
            end
        end
        
        function B = getBMatrix(self, jacobian, der_array, element)
            % Berechnet die Integrationsmatrix in einem Integrationspunkt
            
            % Zunächst werden die Ableitungen nach den x,y,(z)-Koordinaten
            % mithilfe der Jacobimatrix berechnet:
            % dH/dx = J^-1 * dH/dr
            h_derivatives = jacobian\der_array;
            
            number_of_nodes = element.getNumberOfNodes();
            blocks = cell(1,number_of_nodes);
            
            for i=1:number_of_nodes
                h_der = h_derivatives(:,i);
                upper_m = diag(h_der);
                lower_m = self.getPermutations(h_der);
                blocks{i} = vertcat(upper_m, lower_m);
            end

            B = horzcat(blocks{1:end});

        end
        
        function K_C = calcStiffness(self, element)
            K_C = nsNumeric.NumInt.methodIntegrate(@self.constitutiveComponentIntegrator, element.int_points, {});
        end
        
        function F_int = calcLoad(self, element)
            % internal Forces:
            F_int = nsNumeric.NumInt.methodIntegrate(@self.internalForcesIntegrator, element.int_points, struct('element', element) );
            
            % here come the body forces:
            % F_ext = self.calcExternalLoad(element)...
            
            % Transforms the internal-forces-tensortoolbox-tensor back into a matlab-matrix.
            %F = double(F_int); % + F_ext
        end
        
        function BCB = calcBCB(self, element)
            BCB = nsNumeric.NumInt.methodIntegrate(@self.BCBIntegrator, element.int_points,  struct('element', element));
        end
        
    end
    
    methods ( Static )
        function const_comp = constitutiveComponentIntegrator( int_point, ~ )
            B = int_point.B;
            C = int_point.CC;
            dV = int_point.dV;
            const_comp = B' * C * B * dV;
        end
        
        function int_F = internalForcesIntegrator(int_point, parameter)
            int_F = - int_point.B' * int_point.stress * int_point.dV;
        end
        
        function BCB = BCBIntegrator(int_point, parameter)
            BCB = - int_point.B' * int_point.CC * int_point.B * int_point.dV;
        end
            
        function permutations = getPermutations(a)
            if length(a) == 2
                permutations = [a(2), a(1)];
            elseif length(a) == 3
                permutations = [a(2),a(1),0; 0,a(3),a(2); a(3), 0, a(1)];
            end
        end    
        
        function coord = coordinatesToVoigt(A)
            coord = reshape(A',[numel(A),1]);
        end
        
        function M = voigtToMatrix(A)
            if length(A) == 3
                dim = 2;
            elseif length(A) == 6
                dim = 3;
            else
                error('Voigt to Matrix doesn''t work for dim != 2,3')
            end
            
            M = zeros(dim,dim);
            
            if dim == 2
                M(1,1) = A(1);
                M(2,2) = A(2);
                
                M(1,2) = A(3);
                M(2,1) = M(1,2);
            elseif dim == 3
                M(1,1) = A(1);
                M(2,2) = A(2);
                M(3,3) = A(3);
                
                M(1,2) = A(4);
                M(2,3) = A(5);
                M(3,1) = A(6);
                
                M(2,1) = M(1,2);
                M(3,2) = M(2,3);
                M(1,3) = M(3,1);
            else
                error('Voigt to Matrix doesn''t work for dim != 2,3')
            end
        end
    end    
end
