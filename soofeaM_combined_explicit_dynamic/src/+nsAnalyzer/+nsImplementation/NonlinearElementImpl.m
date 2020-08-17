classdef NonlinearElementImpl < nsAnalyzer.nsImplementation.ElementImpl
    methods
        function calcStrainStressInIP(~, element) % make static?
            % In dieser Methode werden für jeden Integrationspunkt des
            % Elements die relevanten Größen berechnet und im Integrationspunkt gespeichert.
            for int_point = element.int_points
%                 disp(['element ',num2str(element.number)])
                
                % Calculate Material Jacobimatrix in int-point
                J = nsAnalyzer.nsJacobian.ElementJacobian(element, int_point, 'material');
                j = nsAnalyzer.nsJacobian.ElementJacobian(element, int_point, 'spatial');
                
                % Hilfsgröße:
                dH = element.type.shape.getDerivativeArray(int_point.natural_int_point.natural_coordinates);
                
                % alternative Berechnung des Deformationsgradienten:
%                 dH_dr = element.type.shape.getDerivativeArray(int_point.natural_int_point.natural_coordinates);
%                 dH_dx = J.getInv()*dH_dr;
%                 du_dx = transpose(dH_dx*element.getDisplacementArray());
                
%                 dim = length(du_dx);
%                 def_grad = eye(dim) + du_dx;
                
                % Erklärung HJ Skript S. 84:
                int_point.HJ = transpose(dH)*transpose(J.getInv());
                int_point.Hj = transpose(dH)*transpose(j.getInv());
                
                % ... und der Rest der relevanten Größen
                int_point.dV = J.getDet();
                int_point.dv = j.getDet();
                
                int_point.F = nsAnalyzer.KinematicTensors.calcDeformationGradient(J, j);
                int_point.C = nsAnalyzer.KinematicTensors.calcRightCauchyGreen(int_point.F);
                int_point.E = nsAnalyzer.KinematicTensors.calcLagrangeStrain(int_point.C);
                
                [S, cc, alpha_new, b_bar_e_new] = element.material.CalcNewStateAndTangentEuler(int_point.F, int_point);
                int_point.S           = S;
                int_point.cc          = cc;
                int_point.alpha_new   = alpha_new;
                int_point.b_bar_e_new = b_bar_e_new;
                
                % calculate cauchy-stress for postprocessing as 3x3 matrix,
                cauchy_stress = 1/det(int_point.F)*int_point.F*int_point.S*transpose(int_point.F);
                int_point.setStress( cauchy_stress );
            end
        end
        
        function K = calcStiffness(self, element)
            % In dieser Funktion wird die Steifigkeitsmatrix eines Elements
            % mittels numerischer Integration ermittelt.
            K_C = nsNumeric.NumInt.methodIntegrate(@self.constitutiveComponentIntegratorEuler, element.int_points, struct('element', element) );
            
            K_S = nsNumeric.NumInt.methodIntegrate(@self.initialStressComponentIntegrator, element.int_points, struct('element', element) );
            
            % Transforms the striffness-matrix-tensortoolbox-tensor back into a matlab-matrix.
            K = double(K_C + K_S);
            
        end
        
        function F = calcLoad(self, element)
            % internal Forces:
            F_int = nsNumeric.NumInt.methodIntegrate(@self.internalForcesIntegrator, element.int_points, struct('element', element) );
            
            % here come the body forces:
            % F_ext = self.calcExternalLoad(element)...
            
            % Transforms the internal-forces-tensortoolbox-tensor back into a matlab-matrix.
            F = double(F_int); % + F_ext
        end
    end
    
    methods ( Static, Access = private )
        function const_comp = constitutiveComponentIntegratorLagrange( int_point, parameter )
            % Weiß der Teufel ob diese Implementierung stimmt...
            % Zumindest die Matrix-Dimensionen müssten stimmen.
            element = parameter.element;
            dimension = element.node_list(1).dimension;
            number_of_nodes = element.getNumberOfNodes();
            dofs_per_element = number_of_nodes * dimension;
            
            CC = element.material.getMaterialElasticityTensor(int_point);
            F = int_point.F;
            HJ = int_point.HJ;
            dV = int_point.dV;
            
            %             temp1 = ttt(HJ,F);
            %             temp2 = ttt(temp1,CC,[2 4],[4 2]);
            %             temp3 = ttt(temp2,F,4,2);
            %             temp4 = ttt(temp3,HJ,3,2);
            
            % Review the indexing process
            %             const_comp_jkil = temp4;
            %             const_comp = permute(const_comp_jkil,[2 3 1 4]);
            
            %             const_comp = double(const_comp);
            
%             temp1 = ttt(HJ,F);
%             temp2 = ttt(temp1,CC,[2 4],[4 2]);
%             temp3 = ttt(temp2,temp1,[3 4],[4 2]);
%             % Review the reshaping process
%             const_comp_array = reshape(temp3,[dofs_per_element dofs_per_element]);
            
            c = zeros(dimension,number_of_nodes,dimension,number_of_nodes);
            for i = 1:dimension
                for j = 1:number_of_nodes
                    for k = 1:dimension
                        for l = 1:number_of_nodes
                            for p = 1:dimension
                                for n = 1:dimension
                                    for m = 1:dimension
                                        for o = 1:dimension
                                            c(i,j,k,l) = c(i,j,k,l) + HJ(j,p)*F(k,n)*CC(m,n,o,p)*F(i,o)*HJ(l,m)*dV;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
         
            const_comp = reshape(c,[dofs_per_element dofs_per_element]);
        end
        
        function const_comp = constitutiveComponentIntegratorEuler( int_point, parameter )
            % Weiß der Teufel ob diese Implementierung stimmt...
            % Zumindest die Matrix-Dimensionen müssten stimmen.
            element = parameter.element;
            dimension = element.node_list(1).dimension;
            number_of_nodes = element.getNumberOfNodes();
            dofs_per_element = number_of_nodes * dimension;
            
            cc = int_point.cc;
            Hj = int_point.Hj;
            dv = int_point.dv;
            
            c = zeros(dimension,number_of_nodes,dimension,number_of_nodes);
            for i = 1:dimension
                for j = 1:number_of_nodes
                    for k = 1:dimension
                        for l = 1:number_of_nodes
                            for m = 1:dimension
                                for n = 1:dimension
                                    c(i,j,k,l) = c(i,j,k,l) + Hj(j,m)*cc(i,m,k,n)*Hj(l,n)*dv;
                                end
                            end
                        end
                    end
                end
            end
         
            
            
            const_comp = reshape(c,[dofs_per_element dofs_per_element]);
        end        
        
        function isc = initialStressComponentIntegrator(int_point, parameter)
            % Weiß der Teufel ob diese Implementierung stimmt...
            % Zumindest die Matrix-Dimensionen müssten stimmen.
            element = parameter.element;
            dimension = element.node_list(1).dimension;
            number_of_nodes = element.getNumberOfNodes();
            dofs_per_element = number_of_nodes * dimension;
            
            HJ = int_point.HJ;
            S = int_point.S;
            dV = int_point.dV;
            delta = eye(dimension);
%             delta = nsNumeric.kronecker(dimension);
            
%             temp1 = ttt(HJ,S,2,2);
%             temp2 = ttt(temp1,nsNumeric.kronecker(dimension));
%             temp3 = ttt(temp2,HJ,2,2);
%             temp4 = temp3*dV;
            
            % Review the indexing process
            
            c = zeros(dimension,number_of_nodes,dimension,number_of_nodes);
            for i = 1:dimension
                for j = 1:number_of_nodes
                    for k = 1:dimension
                        for l = 1:number_of_nodes
                            for n = 1:dimension
                                for m = 1:dimension
                                    c(i,j,k,l) = c(i,j,k,l) + HJ(j,n)*S(m,n)*delta(i,k)*HJ(l,m)*dV;
                                end
                            end
                        end
                    end
                end
            end
            
            % Review the reshaping process
            isc = reshape(c,[dofs_per_element dofs_per_element]);
        end
        
        function int_F = internalForcesIntegrator(int_point, parameter)
            element = parameter.element;
            dimension = element.node_list(1).dimension;
            number_of_nodes = element.getNumberOfNodes();
            dofs_per_element = number_of_nodes * dimension;
            
            HJ = int_point.HJ;
            S = int_point.S;
            F = int_point.F;
            dV = int_point.dV;
            
%             temp1 = ttt(F,S,2,1);
%             temp2 = ttt(temp1,HJ,2,2);
%             
%             int_F = -temp2*dV;
            
            % Indices should already fit.
            
            c = zeros(dimension,number_of_nodes);
            for i = 1:dimension
                for j = 1:number_of_nodes
                    for k = 1:dimension
                        for l = 1:dimension
                            c(i,j) = c(i,j) - F(i,k)*S(k,l)*HJ(j,l)*dV;
                        end
                    end
                end
            end
            
            int_F = reshape(c,[dofs_per_element 1]);
        end
    end
end
