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

classdef LagrangianJacobian < handle
  properties( Access = protected )
    element;
  end
  
  methods
    function self = LagrangianJacobian( element )
      self.element = element;
    end
    
    function determinant = getDet( self, int_point )
      jacobian = self.calc( int_point );
      jac_mat = tenmat( jacobian, 1, 2 );
      determinant = det(jac_mat.data);
    end
    
    function inverse = getInv( self, int_point )
      jacobian = self.calc( int_point );
      jac_mat = tenmat( jacobian, 1, 2 );
      inverse_mat = inv(jac_mat.data);
      inverse = tensor(inverse_mat);
    end
  end
  
  methods( Access = protected )
    function jacobian = calc( self, int_point )
      der_tensor = self.element.element_type.shape.getDerivativeTensor(int_point.reference_coordinates);
      node_point_tensor = self.element.getNodeTensor();
      jacobian = ttt(node_point_tensor,der_tensor,2,1);
    end
  end
end