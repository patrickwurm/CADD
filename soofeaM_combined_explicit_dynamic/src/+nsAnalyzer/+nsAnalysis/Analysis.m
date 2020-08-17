classdef Analysis < handle
    properties
        model
        output_handler
        example
        boundary_cond_matrix %in contrast to standard soofeam, boundary conditions are saved in this matrix, and not in the node DOFs (a lot faster that way)
        number_of_unknowns
        u
        u1
        u_dd
        global_H_matrix_Pad_Atoms
        stiffness_matrix
    end
    
    methods
        function self = Analysis(model, output_handler, example)
            self.model = model;
            self.output_handler = output_handler;
            self.example = example;
            self.number_of_unknowns = self.model.getNumberOfUnknowns();
            self.CalcConstantVariables;
            self.calcStrainStress; %This is done to have the dV needed to compute the mass matrix
            self.stiffness_matrix = self.calcGlobalStiffnessMatrix(zeros(self.number_of_unknowns));
            self.u = zeros(self.number_of_unknowns, 1);
            self.u1 = zeros(self.number_of_unknowns, 1);
            self.u_dd = zeros(self.number_of_unknowns, 1);
        end

        function setModel(self,model)
            self.model = model;
        end
        
        function calcStrainStress(self)
            for element = self.model.element_dict
                element.type.implementation.calcStrainStressInIP(element)
            end
        end
        
        function ComputeHMatrix(self, mxpatm, mxatm)
            self.global_H_matrix_Pad_Atoms = zeros(mxpatm,numel(self.model.node_dict));
            for element = self.model.element_dict
                for pad_atom = element.pad_atoms
                    line_index = pad_atom.number - mxatm;
                    H = element.type.shape.getArray( pad_atom.natural_coordinates );
                    self.global_H_matrix_Pad_Atoms = self.assembleHMatrix(element, H, self.global_H_matrix_Pad_Atoms, line_index);
                end
            end
        end
        
        function CalcConstantVariables(self)
            %Compute everything that is constant over the course of the
            %simulation only ONCE at this point
            for element = self.model.element_dict
                for int_point = element.int_points
                    dimension = length(int_point.material_coordinates);
                    int_point.der_array = element.type.shape.getDerivativeArray(int_point.getNaturalCoordinates());
                    int_point.CC = element.material.getElasticityMatrix(dimension);
                    int_point.jacobian = nsAnalyzer.nsJacobian.Jacobian.staticCalcJacobian(element, int_point, 'material');
                    int_point.dV = det(int_point.jacobian);
                    int_point.B = element.type.implementation.getBMatrix(int_point.jacobian, int_point.der_array, element);
                end
            end
            self.number_of_unknowns = self.model.getNumberOfUnknowns();
        end
        
        function global_stiffness = calcGlobalStiffnessMatrix(self, global_stiffness)
%             disp('|--> Assemble global stiffness')
            
            for element = self.model.element_dict
                % Zunächst wird die Elementsteifigkeitsmatrix in jedem
                % Element berechnet.
                K_elem = element.type.implementation.calcStiffness(element);
                
                % Anschließend wird die Elementsteifigkeitsmatrix in die
                % globale Steifigkeitsmatrix 'hineinassembliert'. Wie diese
                % Methode funktioniert wird im Rahmen der LV nicht
                % behandelt.
                global_stiffness = self.assembleStiffness(element, K_elem, global_stiffness);
            end
        end
        
        function global_load = calcExternalLoadVector(self, global_load)
%NOT NEEDED (ONLY DISPLACEMENT BC IN CADD)
%             disp('|--> Assemble global load')
% 
%             % Body Forces:
%             for element = self.model.element_dict
%                 F_body_element = element.type.implementation.calcBodyForce(element);
%                 global_load = self.assembleLoad(element, F_body_element, global_load);
%             end
% 
% %             % Surface Forces:
%             for boundary = self.model.boundary_dict
%                 for bc = boundary.component_list
%                     F_surface_bc = bc.type.implementation.calcSurfaceLoad( bc );
%                     global_load = self.assembleLoad(bc, F_surface_bc, global_load);
%                 end
%             end
%             
%             % Node Forces:
%             % The node forces dont have to be calculated, they are simply read from each node and written to the according dof in the global load vector.
%             for node = self.model.node_dict
%                 F_node = node.node_force;
%                 if ~all(F_node==0)
%                     for coord_index = 1:self.model.dimension
%                         global_index = (node.number-1) * self.model.dimension + coord_index;
%                         global_load(global_index) = global_load(global_index) + F_node(coord_index);
%                     end
%                 end
%             end
        end
        
        function solveFESystem(self)
            %number of unknowns is misleading, it is the number of all DOF
            %(number of elements * dimension)
            global_load = zeros(self.number_of_unknowns, 1);
            %global_stiffness = zeros(self.number_of_unknowns, self.number_of_unknowns);
            
            % Anschließend werden relevante kinematische und kinetische
            % Größen in allen Integrationspunkten berechnet.
            %self.calcStrainStress();
            
            % Schließlich werden die Elementsteifigkeitsmatrizen
            % berechnet und assembliert.
            %global_stiffness = self.calcGlobalStiffnessMatrix(global_stiffness);
            global_stiffness = self.stiffness_matrix;
            
            % Der globale Kraftvektor wird berechnet und assembliert.
            global_load = self.calcGlobalLoadVector(global_load);
            
            % The Dirichlet boundary conditions are applied to the LSE
            % (Linear System of Equations)
            [global_stiffness, global_load] = self.integrateDirichletBC(global_stiffness, global_load);
                        
            global_stiffness = sparse(global_stiffness);
%             disp('|--> Solve Linear System of Equations')
            % Der Operator '\' löst das GLS
            % global_stiffness*solution_vector = global_load
            %
            % solution_vector are the nodal displacements for a linear
            % analysis and nodal displacement increments for a nonlinear
            % analysis.
            self.u = global_stiffness\global_load;
            
        end
        
        function [global_stiffness, global_load] = integrateDirichletBC(self, global_stiffness, global_load)
            for i=1:length(self.boundary_cond_matrix(:,1))
                        % Zeile im Spaltenvektor der Lösungsvariablen
                        global_col_index = self.boundary_cond_matrix(i,1);
                        bc_value = self.boundary_cond_matrix(i,2);
                        % OLD
%                         for row_counter=1:length(global_stiffness)
%                             if row_counter == global_col_index
%                                 global_load(row_counter) = global_stiffness(row_counter, row_counter)*bc_value;
%                             else
%                                 global_load(row_counter) = global_load(row_counter) - global_stiffness(row_counter, global_col_index)*bc_value;
%                             end
%                         end
                        % NEW
                        global_load = global_load - global_stiffness(:,global_col_index)*bc_value;
                        global_load(global_col_index) = global_stiffness(global_col_index, global_col_index)*bc_value;
                        
                        temp = global_stiffness(global_col_index, global_col_index);
                        global_stiffness(:, global_col_index) = 0;
                        global_stiffness(global_col_index,:) = 0;
                        global_stiffness(global_col_index, global_col_index) = temp;
            end
        end
        
        function paddisp = getPadDisplacements(self)
            displacements = self.u; %+ self.u_dd; %I think this is not the best way of updating the pad
            paddisp = self.global_H_matrix_Pad_Atoms*reshape(displacements',[self.model.dimension,1/2*length(displacements)])';
        end

        % for stress projection
        function calcStrainStressForPostprocessing(self)
            for element = self.model.element_dict
                element.type.implementation.calcStrainStressInIP( element )
                element.type.implementation.calcComparisonStressesInIP( element )
            end
            
            self.performGlobalStressProjection
        end
       
        % for stress projection
        function performGlobalStressProjection(self)
            %stress projection - only for von Mises equivalent stress
            disp('\--> stress projection - only for von Mises equivalent stress')
            if( self.model.dimension == 2 || self.model.dimension == 3 )
                stress_projection_vector = zeros( self.model.getNumberOfNodes,1 );
                mass_matrix = zeros( self.model.getNumberOfNodes, self.model.getNumberOfNodes );
            end
            for element = self.model.element_dict
                stress_projection_vector = self.assembleLoad(element, element.type.implementation.calcProjectionVector(element),...
                    stress_projection_vector);
                mass_matrix = self.assembleStiffness(element, element.type.implementation.calcMassMatrix( element ),...
                    mass_matrix);
            end
            nodal_stress_vector = mass_matrix\stress_projection_vector;
            for node = self.model.node_dict
                node.setStressVonMises(nodal_stress_vector(node.number));
            end
        end        
    end
    
    methods( Static )
        function global_stiffness = assembleStiffness(node_container, local_stiffness, global_stiffness)
            
            % dofs_per_node is dimension of stiffness matrix divided by number of nodes per element.
            % I have not figured out why, yet. Maybe number of entries
            % means numbers of dofs per node, so 2 for 2D and 3 for 3D
            % original line: nr_of_entries = local_stiffness.shape(1)/node_container.type.shape.getNumberOfNodes;
            
            dofs_per_element = length(local_stiffness);
            
            dofs_per_node = dofs_per_element/node_container.type.shape.getNumberOfNodes;
            
            for local_row_index = 1:length(local_stiffness)
                for local_col_index = 1:length(local_stiffness)
                    % In the assembly process, the task is to find the
                    % global indices of the global stiffness matrix which
                    % belong to the local indices of the local indices of
                    % the current element
                    
                    % There are dofs_per_node columns for each node.
                    % example 2D: u1,v1,u2,v2,...
                    % -> indices 1 & 2 must be mapped to 1
                    % -> indices 3 & 4 must be mapped to 2, etc.
                    % hence the ceil.
                    local_node_index_row = ceil(local_row_index/dofs_per_node);
                    local_node_index_col = ceil(local_col_index/dofs_per_node);
                    
                    % get the current nodes
                    row_node = node_container.node_list(local_node_index_row);
                    col_node = node_container.node_list(local_node_index_col);
                    
                    % get the global node number
                    global_node_index_row = row_node.number;
                    global_node_index_col = col_node.number;
                    
                    % retrieve the position in the global stiffness matrix.
                    % example 2D: node 8 with dofs u8 & v8 must deliver
                    % global indices 7*2+1 = 15 & 7*2+2 = 16
                    % the first summend is called start_position, the
                    % second one is called offset and can be obtained via a
                    % modulu operation.
                    row_start_position = (global_node_index_row - 1)*dofs_per_node;
                    col_start_position = (global_node_index_col - 1)*dofs_per_node;
                    
                    row_offset = mod(local_row_index,dofs_per_node);
                    col_offset = mod(local_col_index,dofs_per_node);
                    
                    
                    % at the last index of each node, the mod operator
                    % returns 0, but it should of course return
                    % dofs_per_node
                    if row_offset == 0
                        row_offset = dofs_per_node;
                    end
                    if col_offset == 0
                        col_offset = dofs_per_node;
                    end
                    
                    % so the global indices are finally obtained
                    global_row_index = row_start_position + row_offset;
                    global_col_index = col_start_position + col_offset;
                    
                    global_stiffness(global_row_index,global_col_index) = global_stiffness(global_row_index,global_col_index) + local_stiffness(local_row_index,local_col_index);
                end
            end
        end
        
        function global_load = assembleLoad(node_container, local_load, global_load)
            
            % analogous to assembleStiffness(...) method above
            dofs_per_node = length(local_load)/node_container.type.shape.getNumberOfNodes;
            for local_row_index = 1:length(local_load)
                
                local_node_index_row = ceil(local_row_index/dofs_per_node);
                row_node = node_container.node_list(local_node_index_row);
                
                global_node_index_row = row_node.number;
                
                row_start_position = (global_node_index_row - 1)*dofs_per_node;
                row_offset = mod(local_row_index,dofs_per_node);
                if row_offset == 0
                    row_offset = dofs_per_node;
                end
                global_row_index = row_start_position + row_offset;
                
                global_load(global_row_index) = global_load(global_row_index) + local_load(local_row_index);
            end
        end
        
        function global_H = assembleHMatrix(node_container, local_H, global_H, line_index)
            %This assembles the H matrix needed to update the Pad atom
            %displacements via the global node displacements
            
            for i=1:numel(node_container.node_list)
                col_index = node_container.node_list(i).number;
                global_H(line_index, col_index) = local_H(i);
            end
            
        end
        
    end

    methods ( Abstract )
        run(self)
    end
        
    methods ( Abstract, Access = protected )
        updateDOF(self, solution_vector)
        calcGlobalLoadVector(self, global_load)
    end
end

