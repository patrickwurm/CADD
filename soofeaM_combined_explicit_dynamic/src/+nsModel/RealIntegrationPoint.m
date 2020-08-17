classdef RealIntegrationPoint < handle
    % Each nsModel.NodeContainer gets its own set of
    % soofeaM.nsModel.RealIntegrationPoint objects
    
    % set-functions need to be implemented!
    properties % ( SetAccess = protected )
        % Each Integration point has an associated natural integration
        % point on the reference cell.
        natural_int_point
        
        % For a lagrange formulation, it suffices to know the material
        % coordinates of each integration point.
        material_coordinates = []
        
        % general quantities used for linear and nonlinear formulations
        stress           % Cauchy-Spannung im Integrationspunkt
        stress_von_mises % Von-Mises-Vergleichsspannung
        body_force       % Volumskräfte
        surface_force    % Oberflächenkräfte
        dV               % enthält die material Jacobi-Determinante
        dv               % enthält die spatial Jacobi-Determinante
    
        % quantities used in a linear formulation
        B      % B-Matrix im Integrationspunkt in Voigt-Notation
        strain % Verzerrung im Integrationspunkt
        CC     % Elastizitätsmatrix in Voigt-Notation
        
        % quantities used in a nonlinear formulation
        HJ    % see p.84 im Hammer-Skriptum
        Hj    % Äquivalent dazu mit Euler-Jacobi-Determinante
        F     % deformation gradient
        C     % right cauchy-green deformation tensor
        E     % green-lagrange strain tensor
        S     % second piola-kirchhoff stress tensor
        
        cc    % spatial material tangent
        
        % quantities for a hyperplastic formulation
        F_old
        b_bar_e_old
        alpha_old
        
        b_bar_e_new
        alpha_new
        
        der_array
        jacobian
    end
    
    methods
        function self = RealIntegrationPoint( natural_int_point, material_coordinates )
            self.natural_int_point = natural_int_point;
            self.material_coordinates = material_coordinates;
            
            dim = length(material_coordinates);
            
            self.F_old       = eye(dim);
            self.b_bar_e_old = eye(dim);
            self.alpha_old   = 0;
        end
        
        function weight = getWeight(self)
            weight = self.natural_int_point.weight;
        end
        
        function natural_coords = getNaturalCoordinates(self)
            natural_coords = self.natural_int_point.natural_coordinates;
        end
        
%== A bunch of set-functions missing. maybe make structs for properties ==%
        function setB(self, B)
            self.B = B;
        end
        
        function setStrain(self, strain)
            self.strain = strain;
        end
        
        function setCC(self, CC)
            self.CC = CC;
        end        
        
        function setStress(self,stress)
            self.stress = stress;
        end
        
        function setSurfaceForce(self, surface_force)
            self.surface_force = surface_force;
        end
        
        function setBodyForce(self, body_force)
            self.body_force = body_force;
        end
        
        function setStressVonMises(self, svm)
            self.stress_von_mises = svm;
        end
    end
end
