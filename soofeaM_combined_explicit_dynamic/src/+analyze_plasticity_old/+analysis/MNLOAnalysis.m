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

classdef MNLOAnalysis < analyze.analysis.Analysis
  methods
    function self = MNLOAnalysis( model , model_progression )
      self@analyze.analysis.Analysis( model , model_progression );
    end
    
    function run(self, data_set)
      disp('> Material Nonlinearity Only Analysis (MNLO)');
      
      convergence_steps = data_set.number_of_convergence_steps;
      iterations = data_set.number_of_iteration_steps_global;
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %VTK OUTPUT%
      global global_data;
      sim_name = global_data.sim_name;
      filestr = strcat('examples/',sim_name,'.vtk.*');
      delete(filestr)
      self.VTKout(0)
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      for conv=1:convergence_steps
          
          disp(['Step: ', num2str(conv),'/',num2str(convergence_steps), ...
              '-------------------------------------------------------------',...
              'Step: ', num2str(conv),'/',num2str(convergence_steps) ]);
      
          for iter=1:iterations
              %integrate boundary conditions
              self.bc_handler.integrateBC(iter, conv, convergence_steps);
              
              % find Delta of DOF
              load_vector = [];
              global_stiffness = [];
              [global_stiffness,load_vector] = self.assemble_k_fint__(global_stiffness,load_vector,iter);
              [global_stiffness,load_vector] = self.assemble_fext__(global_stiffness,load_vector);
              
              solution = global_stiffness \ load_vector;
              
              % update DOF
              self.updateDOF__( solution );
              
              % update Eulerian coordinates
              for node_counter = 1:self.model.nodeAmount()
                  self.model.getNode(node_counter).updateEulerianCoordUnconstrained();
                  if (iter==1)
                      self.model.getNode(node_counter).updateEulerianCoordConstrained();
                  end
              end
              %  for i=1:self.model.nodeAmount()
              %      self.model.getNode( i ).dispIterationDisplacement();
              %  end
              
              % check abort criterion
              unbalanced_energy = self.unbalanced_energy__(solution, load_vector);
              disp(['--- <a href="">Unbalanced Energy</a> (iter=', num2str(iter), ...
                  '): ', num2str(unbalanced_energy)]);
              
              if ( unbalanced_energy < data_set.tolerance_global )
                  break;
              end
              
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if (2<1)
                  disp('Visualisation: ');
                  visualisation = visualise.Visualisation( self.model );
                  visualisation.visualise( conv );
              end
              %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              
          end
          
          for el_counter = 1:self.model.elementAmount()
              for i=1:numel(self.model.getElement(el_counter).integration_points);
                self.model.getElement(el_counter).integration_points(i).SetStateOld( self.model.getElement(el_counter).integration_points(i).GetStateNew);
                sn = self.model.getElement(el_counter).integration_points(i).GetStressNew;
                self.model.getElement(el_counter).integration_points(i).SetStressOld( sn );
                
                delta = tensor(math.kronecker([3 3]));
                
                devs = sn - 1/3 * ttt(sn,delta,[1 2], [1 2]) * delta;
                svM = sqrt(3/2)*norm(devs);
                self.model.getElement(el_counter).integration_points(i).SetsvM( svM );
                
                en = self.model.getElement(el_counter).integration_points(i).GetStrain;
                deve = en - 1/3 * ttt(en,delta,[1 2], [1 2]) * delta;
                evM = sqrt(2/3)*norm(deve);
                self.model.getElement(el_counter).integration_points(i).SetevM( evM );
               end
          end
          
          for node_counter = 1:self.model.nodeAmount()
              self.model.getNode(node_counter).updateLastConvergedCoordUnconstrained();
              self.model.getNode(node_counter).updateLastConvergedCoordConstrained();
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Projection of Intergration Point Variables to Nodes
          
          %Clear nodal stresses, strains and internal variables from
          %previous increment
          for node_counter = 1:self.model.nodeAmount()
            self.model.getNode(node_counter).resetStressStrainInt(); 
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %STRESSES
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %loop over all elements
          for el_counter = 1:self.model.elementAmount()
              %get "hat" vector
              
              a = zeros(numel(self.model.getElement(el_counter).integration_points),1);
              
              %loop over all nodes and get the ansatz function vector
              for i=1:numel(self.model.getElement(el_counter).node_array);
                  myStress = zeros(3,3);
                  H = self.model.getElement(el_counter).element_type.shape.getTensorExtrapolation(self.model.getElement(el_counter).element_type.shape.getReferenceCoordinates(i));
                  for str_i = 1:3
                      for str_j = str_i:3
                          for j=1:numel(self.model.getElement(el_counter).integration_points)
                              a(j,1) = self.model.getElement(el_counter).integration_points(j).stress_new(str_i,str_j);
                          end
                          %hat vector needs to be resorted, because H is sorted in nodal ordering 
                          %while hat vector is sorted in int_point
                          %ordering!
                          b = zeros(numel(self.model.getElement(el_counter).integration_points),1);
                          b(1,1) = a(1,1);
                          b(2,1) = a(3,1);
                          b(3,1) = a(7,1);
                          b(4,1) = a(5,1);
                          b(5,1) = a(2,1);
                          b(6,1) = a(4,1);
                          b(7,1) = a(8,1);
                          b(8,1) = a(6,1);
                          myStress(str_i,str_j) =  dot(H,b);
                          if str_i~=str_j
                            myStress(str_j,str_i) = myStress(str_i,str_j);
                          end
                      end
                  end
                  if conv==1
                    self.model.getElement(el_counter).node_array(i).addNeighbor( );
                  end
                  self.model.getElement(el_counter).node_array(i).setStress ( myStress );
              end  
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %STRAINS
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %loop over all elements
          for el_counter = 1:self.model.elementAmount()
              %get "hat" vector
              
              a = zeros(numel(self.model.getElement(el_counter).integration_points),1);
              
              %loop over all nodes and get the ansatz function vector
              for i=1:numel(self.model.getElement(el_counter).node_array);
                  myStrain = zeros(3,3);
                  H = self.model.getElement(el_counter).element_type.shape.getTensorExtrapolation(self.model.getElement(el_counter).element_type.shape.getReferenceCoordinates(i));
                  for str_i = 1:3
                      for str_j = str_i:3
                          for j=1:numel(self.model.getElement(el_counter).integration_points)
                              a(j,1) = self.model.getElement(el_counter).integration_points(j).strain(str_i,str_j);
                          end
                          %hat vector needs to be resorted, because H is sorted in nodal ordering
                          %while hat vector is sorted in int_point
                          %ordering!
                          b = zeros(numel(self.model.getElement(el_counter).integration_points),1);
                          b(1,1) = a(1,1);
                          b(2,1) = a(3,1);
                          b(3,1) = a(7,1);
                          b(4,1) = a(5,1);
                          b(5,1) = a(2,1);
                          b(6,1) = a(4,1);
                          b(7,1) = a(8,1);
                          b(8,1) = a(6,1);
                          myStrain(str_i,str_j) =  dot(H,b);
                          if str_i~=str_j
                              myStrain(str_j,str_i) = myStrain(str_i,str_j);
                          end
                      end
                  end
                  self.model.getElement(el_counter).node_array(i).setStrain ( myStrain );
              end  
          end 
          
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           PLASTIC STRAINS (not needed, rather go for equivalent plastic strain)
%           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           loop over all elements
%           for el_counter = 1:self.model.elementAmount()
%               get "hat" vector
%               
%               a = zeros(numel(self.model.getElement(el_counter).integration_points),1);
%               
%               loop over all nodes and get the ansatz function vector
%               for i=1:numel(self.model.getElement(el_counter).node_array);
%                   myStrain = zeros(3,3);
%                   H = self.model.getElement(el_counter).element_type.shape.getTensorExtrapolation(self.model.getElement(el_counter).element_type.shape.getReferenceCoordinates(i));
%                   for str_i = 1:3
%                       for str_j = str_i:3
%                           for j=1:numel(self.model.getElement(el_counter).integration_points)
%                               a(j,1) = self.model.getElement(el_counter).integration_points(j).state_new{2}(str_i,str_j);
%                           end
%                           hat vector needs to be resorted, because H is sorted in nodal ordering
%                           while hat vector is sorted in int_point
%                           ordering!
%                           b = zeros(numel(self.model.getElement(el_counter).integration_points),1);
%                           b(1,1) = a(1,1);
%                           b(2,1) = a(3,1);
%                           b(3,1) = a(7,1);
%                           b(4,1) = a(5,1);
%                           b(5,1) = a(2,1);
%                           b(6,1) = a(4,1);
%                           b(7,1) = a(8,1);
%                           b(8,1) = a(6,1);
%                           myStrain(str_i,str_j) =  dot(H,b);
%                           if str_i~=str_j
%                               myStrain(str_j,str_i) = myStrain(str_i,str_j);
%                           end
%                       end
%                   end
%                   self.model.getElement(el_counter).node_array(i).setPlasticStrain ( myStrain );
%               end  
%           end 

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %OTHER STATE VARIABLES
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %loop over all elements
          for el_counter = 1:self.model.elementAmount()
              %get "hat" vector
              
              a = zeros(numel(self.model.getElement(el_counter).integration_points),1);
              
              %loop over all nodes and get the ansatz function vector
              for i=1:numel(self.model.getElement(el_counter).node_array);
                  myIntVar = 0;
                  H = self.model.getElement(el_counter).element_type.shape.getTensorExtrapolation(self.model.getElement(el_counter).element_type.shape.getReferenceCoordinates(i));
                  for j=1:numel(self.model.getElement(el_counter).integration_points)
                      a(j,1) = self.model.getElement(el_counter).integration_points(j).state_new{1,1};
                  end
                  %hat vector needs to be resorted, because H is sorted in nodal ordering
                  %while hat vector is sorted in int_point
                  %ordering!
                  b = zeros(numel(self.model.getElement(el_counter).integration_points),1);
                  b(1,1) = a(1,1);
                  b(2,1) = a(3,1);
                  b(3,1) = a(7,1);
                  b(4,1) = a(5,1);
                  b(5,1) = a(2,1);
                  b(6,1) = a(4,1);
                  b(7,1) = a(8,1);
                  b(8,1) = a(6,1);
                  myIntVar =  dot(H,b);
                  self.model.getElement(el_counter).node_array(i).setAlpha(myIntVar);
              end
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %VON MISES STRESS
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %loop over all elements
          for el_counter = 1:self.model.elementAmount()
              %get "hat" vector
              
              a = zeros(numel(self.model.getElement(el_counter).integration_points),1);
              
              %loop over all nodes and get the ansatz function vector
              for i=1:numel(self.model.getElement(el_counter).node_array);
                  svM = 0;
                  H = self.model.getElement(el_counter).element_type.shape.getTensorExtrapolation(self.model.getElement(el_counter).element_type.shape.getReferenceCoordinates(i));
                  for j=1:numel(self.model.getElement(el_counter).integration_points)
                      a(j,1) = self.model.getElement(el_counter).integration_points(j).svM;
                  end
                  %hat vector needs to be resorted, because H is sorted in nodal ordering
                  %while hat vector is sorted in int_point
                  %ordering!
                  b = zeros(numel(self.model.getElement(el_counter).integration_points),1);
                  b(1,1) = a(1,1);
                  b(2,1) = a(3,1);
                  b(3,1) = a(7,1);
                  b(4,1) = a(5,1);
                  b(5,1) = a(2,1);
                  b(6,1) = a(4,1);
                  b(7,1) = a(8,1);
                  b(8,1) = a(6,1);
                  svM =  dot(H,b);
                  self.model.getElement(el_counter).node_array(i).setsvM ( svM );
              end
          end
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %VON MISES STRAIN (equivalent plastic strain)
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %loop over all elements
          for el_counter = 1:self.model.elementAmount()
              %get "hat" vector
              
              a = zeros(numel(self.model.getElement(el_counter).integration_points),1);
              
              %loop over all nodes and get the ansatz function vector
              for i=1:numel(self.model.getElement(el_counter).node_array);
                  evM = 0;
                  H = self.model.getElement(el_counter).element_type.shape.getTensorExtrapolation(self.model.getElement(el_counter).element_type.shape.getReferenceCoordinates(i));
                  for j=1:numel(self.model.getElement(el_counter).integration_points)
                      a(j,1) = self.model.getElement(el_counter).integration_points(j).evM;
                  end
                  %hat vector needs to be resorted, because H is sorted in nodal ordering
                  %while hat vector is sorted in int_point
                  %ordering!
                  b = zeros(numel(self.model.getElement(el_counter).integration_points),1);
                  b(1,1) = a(1,1);
                  b(2,1) = a(3,1);
                  b(3,1) = a(7,1);
                  b(4,1) = a(5,1);
                  b(5,1) = a(2,1);
                  b(6,1) = a(4,1);
                  b(7,1) = a(8,1);
                  b(8,1) = a(6,1);
                  evM =  dot(H,b);
                  self.model.getElement(el_counter).node_array(i).setevM ( evM );
              end
          end
          
          %AVERAGE ALL NODAL QUANTITIES (DIVIDE SUM BY AMOUNT OF NEIGHBORING ELEMENTS)
          for node_counter = 1:self.model.nodeAmount()
            node = self.model.getNode(node_counter);
            node.stress = node.stress/node.nneighbor_elements;
            node.strain = node.strain/node.nneighbor_elements;
            node.svM = node.svM/node.nneighbor_elements;
            node.evM = node.evM/node.nneighbor_elements;
            node.internal_variables{1,1} = node.internal_variables{1,1}/node.nneighbor_elements;
            node.internal_variables{1,2} = node.internal_variables{1,2}/node.nneighbor_elements;
          end
          
          %Plot stresses just to check, can be deleted
%           for el_counter = 1:self.model.elementAmount()
%               disp(['Element #',num2str(el_counter)])
%               for i=1:numel(self.model.getElement(el_counter).node_array);
%                   disp(['NODAL STRESSES, NODE #',num2str(i),'Coords:',num2str(self.model.getElement(el_counter).node_array(i).lagrangian_coords)])
%                   self.model.getElement(el_counter).node_array(i).stress
%                   disp(['INT_POINT STRESSES, IP #',num2str(i)])
%                   self.model.getElement(el_counter).integration_points(i).stress_new
%                   disp('##############################################################')
%               end
%           end

          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %VTK OUTPUT%
          
          self.VTKout(conv)
          

      end
    end    
    
    function VTKout(self, conv)
        global global_data;
        filestr = strcat('examples/',global_data.sim_name,'.vtk.');
        outputfile = strcat(filestr,num2str(conv));
        dlmwrite(outputfile,  '# vtk DataFile Version 3.0', 'delimiter', '', '-append');
        dlmwrite(outputfile,  'VTK Output', 'delimiter', '', '-append');
        dlmwrite(outputfile,  'ASCII', 'delimiter', '', '-append');
        dlmwrite(outputfile,  'DATASET UNSTRUCTURED_GRID', 'delimiter', '', '-append');
        dlmwrite(outputfile,  strjoin({'POINTS',num2str(self.model.nodeAmount()),'float'}), 'delimiter', '', '-append');
        for node_counter = 1:self.model.nodeAmount()
            out = self.model.getNode(node_counter).lagrangian_coords;
            dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
        end
        %dlmwrite(outputfile,  ' ', 'delimiter', '', '-append');
        dlmwrite(outputfile,  strjoin({'CELLS',num2str(self.model.elementAmount()),num2str(9*self.model.elementAmount())}), 'delimiter', '', '-append');
        for el_counter = 1:self.model.elementAmount()
            %Reorder Node Indexing to be consistent with VTK
            a = [self.model.getElement(el_counter).node_array.number];
            a = a-1;
            a = flipud(reshape(a,2,[]));
            b = a(1:2,1:2);
            a(1:2,1:2)=a(1:2,3:4);
            a(1:2,3:4)=b;
            a = reshape(a,1,[]);
            out = [8 a];
            dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
        end
        dlmwrite(outputfile,  ' ', 'delimiter', '', '-append');
        dlmwrite(outputfile,  strjoin({'CELL_TYPES',num2str(self.model.elementAmount())}), 'delimiter', '', '-append');
        for el_counter = 1:self.model.elementAmount()
            dlmwrite(outputfile,  '12', 'delimiter', '', '-append');
        end
        dlmwrite(outputfile,  ' ', 'delimiter', '', '-append');
        dlmwrite(outputfile,  strjoin({'POINT_DATA',num2str(self.model.nodeAmount())}), 'delimiter', '', '-append');
        dlmwrite(outputfile,  'VECTORS displacement double', 'delimiter', '', '-append');
        for node_counter = 1:self.model.nodeAmount()
            out = self.model.getNode(node_counter).eulerian_coords-self.model.getNode(node_counter).lagrangian_coords;
            dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
        end
        dlmwrite(outputfile,  'TENSORS stresses double', 'delimiter', '', '-append');
        for node_counter = 1:self.model.nodeAmount()
            for j=1:3
                out = [ self.model.getNode(node_counter).stress(j,1), self.model.getNode(node_counter).stress(j,2), self.model.getNode(node_counter).stress(j,3) ];
                dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
            end
        end
        dlmwrite(outputfile,  'TENSORS strains double', 'delimiter', '', '-append');
        for node_counter = 1:self.model.nodeAmount()
            for j=1:3
                out = [ self.model.getNode(node_counter).strain(j,1), self.model.getNode(node_counter).strain(j,2), self.model.getNode(node_counter).strain(j,3) ];
                dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
            end
        end
        dlmwrite(outputfile,  'SCALARS alpha double 1', 'delimiter', '', '-append');
        dlmwrite(outputfile,  'LOOKUP_TABLE default', 'delimiter', '', '-append');
        for node_counter = 1:self.model.nodeAmount()
            out = [ self.model.getNode(node_counter).internal_variables{1,1}];
            dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
        end
        dlmwrite(outputfile,  'SCALARS svM double 1', 'delimiter', '', '-append');
        dlmwrite(outputfile,  'LOOKUP_TABLE default', 'delimiter', '', '-append');
        for node_counter = 1:self.model.nodeAmount()
            out = [ self.model.getNode(node_counter).svM];
            dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
        end
        dlmwrite(outputfile,  'SCALARS evM double 1', 'delimiter', '', '-append');
        dlmwrite(outputfile,  'LOOKUP_TABLE default', 'delimiter', '', '-append');
        for node_counter = 1:self.model.nodeAmount()
            out = [ self.model.getNode(node_counter).evM];
            dlmwrite(outputfile, out, 'delimiter', ' ', 'precision', 15, '-append');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
  end
end