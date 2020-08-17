classdef LinearHexagonalMaterial < nsModel.nsMaterial.Material

    properties
        c11
        c12
        c66
        rho
        nu
        mu
    end
    
    methods
        function self = LinearHexagonalMaterial(number, c11, c12, rho)        
            self@nsModel.nsMaterial.Material(number);
            if nargin < 4
                self.rho = 0;
            else
                self.rho = rho;
            end
            self.c11 = c11;
            self.c12 = c12;
            self.c66 = 1/2*(self.c11-self.c12);
            [self.mu, self.nu] = self.computeVoigtAverageElasticModuli(self.c11,self.c12,self.c66);
        end
        
        function C = getElasticityMatrix(self,dimension)
            if dimension == 2
                C = [self.c11, self.c12, 0; self.c12, self.c11, 0; 0, 0, self.c66];
            elseif dimension == 3
                error('3 dimensions not supported')
            end
        end
        
        function [mu, nu] = computeVoigtAverageElasticModuli(self, c11, c12, c66)
            %see Hirth and Lothe Page 430 Cubic Crystals (I hope this is right)
            H = 2*c66 + c12 - c11;
            mu = c66 - 1/5*H; % = G = shear modulus
            lambda = c12 - 1/5*H;
            nu = lambda / (2*lambda + mu); % See Wikipedia Lame Parameters
            % given: lambda, G; wanted: nu
        end
        
        

    end
    
end

