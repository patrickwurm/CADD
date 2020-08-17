classdef Explicit_Dynamic_Analysis < handle
    properties
        model
        output_handler
        damp_matrix
        mass_matrix % Should the mass matrix be another property of Dynamic_Analysis?
        u % node displacements; They are stored at the node level in the standard version of soofeam,
        % however, storing them in a global vector is a lot faster; when an
        % output is requested, the current displacements are written to
        % the node objects
        u_dd %nodal displacements in to the infinite continuum dislocation region, u_total = u + u_dd
        u1 % node velocities; They should actually also be stored at the node level 
        %to comply with standard soofeam; but that is slow
        u2 % node accelerations; see node velocities
        lumped_mass_matrix %scalar value to indicate the use of a lumped mass matrix; is either 1 or 0
        global_BCB_matrix %only calculated once at the beginning of the simulation 
        % (this quantity plays a big role in the vectorization of the material routine/internal node force calculation)
        number_of_unknowns %isnt actually the number of unknowns, but the number of all degrees of freedom
        boundary_cond_matrix %in contrast to standard soofeam, boundary conditions are saved in this matrix, and not in the node DOFs (a lot faster that way)
        global_H_matrix_Pad_Atoms
        example
        DD_nodal_coordinates
        include_external_BCs %this is set to 0 for the hybrid continuum, to set the internal BCs to zero
    end
    
    methods
        function self = Explicit_Dynamic_Analysis(model, output_handler, example, include_external_BCs)
            self.model = model;
            self.output_handler = output_handler;
            number_of_unknowns = self.model.getNumberOfUnknowns();
            self.CalcConstantVariables;
            self.calcStrainStress; %This is done to have the dV needed to compute the mass matrix
            self.mass_matrix = self.ComputeMassMatrix(zeros(number_of_unknowns));
            self.damp_matrix = self.mass_matrix/self.model.material_dict{1,1}.rho * 0; 
            self.global_BCB_matrix = self.ComputeBCBMatrix(zeros(number_of_unknowns));
            self.u = zeros(number_of_unknowns, 1);
            self.u_dd = zeros(number_of_unknowns, 1);
            self.u1 = zeros(number_of_unknowns, 1);
            self.u2 = zeros(number_of_unknowns, 1);
            self.example = example;
            self.DD_nodal_coordinates = [];
            self.number_of_unknowns = number_of_unknowns;
            self.include_external_BCs = include_external_BCs;
        end

        function setModel(self,model)
            self.model = model;
        end
        
        function calcStrainStress(self)
            for element = self.model.element_dict
                element.type.implementation.calcStrainStressInIP(element)
            end
        end
        
        function mass_matrix = ComputeMassMatrix(self, mass_matrix)

            for element = self.model.element_dict
                % Calculate element mass matrix
                M_elem = 1/2 * self.model.material_dict{1,1}.rho * element.type.implementation.calcMassMatrix(element,1);
                
                % Lumped element mass matrix via row sum method:
                M_elem_lumped = zeros(size(M_elem));
                
                
                for i=1:length(M_elem_lumped(1,:))
                    M_elem_lumped(i,i) = sum(M_elem(i,:));
                end
                
                % Assemble lumped mass matrix into global one
                mass_matrix = self.assembleMassMatrix(element, M_elem_lumped, mass_matrix);
            end
            self.lumped_mass_matrix = 1; %1; we use a lumped mass matrix
        end
        
        function global_BCB_matrix = ComputeBCBMatrix(self, global_BCB_matrix)
            
            for element = self.model.element_dict
                BCB = element.type.implementation.calcBCB(element);
                global_BCB_matrix = self.assembleMassMatrix(element, BCB, global_BCB_matrix);
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
        
        %NOT NEEDED (ONLY DISPLACEMENT BC IN CADD)
        function global_load = calcExternalLoadVector(self, global_load)

            % Node Forces:
            % The node forces dont have to be calculated, they are simply read from each node and written to the according dof in the global load vector.
            for node = self.model.node_dict
                F_node = node.node_force;
                if ~all(F_node==0)
                    for coord_index = 1:self.model.dimension
                        global_index = (node.number-1) * self.model.dimension + coord_index;
                        global_load(global_index) = global_load(global_index) + F_node(coord_index);
                    end
                end
            end
        end
        
        function global_load = calcGlobalLoadVector_vectorized(self, global_load)
            global_load = global_load + self.global_BCB_matrix * self.u;
        end
        
        function solveFESystem(self)
            
            %number of unknowns is misleading, it is the number of all DOF
            %(number of elements * dimension)
            global_load = zeros(self.number_of_unknowns, 1);

            %CALCULATE NEW DISPLACEMENTS
            self.u = self.u + self.u1 * self.model.time_step + 1/2 * self.u2 * self.model.time_step^2;
            
            %INCORPORATE KNOWN DISPLACEMENTS
            self.u = self.integrateDisplacementorOrVelocityBC(self.u, 'disp');
            
            %estimate of v
            %u1_estimate = self.u1 + self.u2 * self.model.time_step; %there is no 1/2 in this estimate
            %u1_estimate = self.integrateDisplacementorOrVelocityBC(u1_estimate, 'vel'); %must this be done?
            
            %LATER
            % Incorporate Atomic Displacements
            %self.incorporateAtomicDisplacements(interface_atom_displacements);
            
            %Compute the load vector in a vectorzied manner
            global_load = self.calcGlobalLoadVector_vectorized(global_load);

            % Incorporate Atomic Accelerations into the LSOE
            [adapted_mass_matrix, global_load] = self.integrateAccelerationBC(self.mass_matrix, global_load);
            
            adapted_mass_matrix = sparse(adapted_mass_matrix);

            % note that the solution vector here consists of accelerations
            % (in contrast to displacements in standard soofeam)
            solution_vector = adapted_mass_matrix\global_load;

            u1_estimate = self.u1 + self.model.time_step*self.u2;
            u1_estimate = self.integrateDisplacementorOrVelocityBC(u1_estimate, 'vel');
            
            u1_hat = self.u1 + 1/2 * self.model.time_step * ( self.u2 + solution_vector - adapted_mass_matrix\(self.damp_matrix*u1_estimate));
            u1_hat = self.integrateDisplacementorOrVelocityBC(u1_hat, 'vel');
            
            a_hat = solution_vector - adapted_mass_matrix\(self.damp_matrix*u1_hat);
            
            self.u1 = self.u1 + 1/2 * self.model.time_step * (self.u2 + a_hat);
            self.u1 = self.integrateDisplacementorOrVelocityBC(self.u1, 'vel');
            
            self.u2 = solution_vector;
        end


        function [global_mass, global_load] = integrateAccelerationBC(self, global_mass, global_load)
            % This has the exact same function as the
            % "integrateDirichletBC" function in standard soofeam
            % (although it becomes slightly simpler due to the fact that the lumped mass matrix is a diagonal matrix)
            % This was optimized quite a bit, we do not go through all
            % nodes anymore to find the constrained ones, but find the
            % constrained ones ONCE at the beginning and safe the current
            % boundary conditions inside a simple array (together with global indices)
            
            for i=1:length(self.boundary_cond_matrix(:,1))
                global_load(self.boundary_cond_matrix(i,1)) = global_mass(self.boundary_cond_matrix(i,1),self.boundary_cond_matrix(i,1))*self.boundary_cond_matrix(i,4);
            end
            
        end
        
        function disp_vel_vector = integrateDisplacementorOrVelocityBC(self, disp_vel_vector, type)
            % See integrateAccelerationBC
            if strcmp(type,'disp')
                index = 2;
            elseif strcmp(type,'vel')
                index = 3;
            end
            for i=1:length(self.boundary_cond_matrix(:,1))
                disp_vel_vector(self.boundary_cond_matrix(i,1)) = self.boundary_cond_matrix(i,index);
            end
        end
        
        function paddisp = getPadDisplacements(self)
            displacements = self.u; %+ self.u_dd; %I think this is not the best way of updating the pad
            paddisp = self.global_H_matrix_Pad_Atoms*reshape(displacements',[self.model.dimension,1/2*length(displacements)])';
        end
        
        %NOT NEEDED (We don't need no stress projection YET?)
%         % for stress projection
%         function calcStrainStressForPostprocessing(self)
%             for element = self.model.element_dict
%                 element.type.implementation.calcStrainStressInIP( element )
%                 element.type.implementation.calcComparisonStressesInIP( element )
%             end
%             
%             self.performGlobalStressProjection
%         end

        %NOT NEEDED (We don't need no stress projection YET?)
%         % for stress projection
%         function performGlobalStressProjection(self)
%             %stress projection - only for von Mises equivalent stress
%             disp('\--> stress projection - only for von Mises equivalent stress')
%             if( self.model.dimension == 2 || self.model.dimension == 3 )
%                 stress_projection_vector = zeros( self.model.getNumberOfNodes,1 );
%                 mass_matrix = zeros( self.model.getNumberOfNodes, self.model.getNumberOfNodes );
%             end
%             for element = self.model.element_dict
%                 stress_projection_vector = self.assembleLoad(element, element.type.implementation.calcProjectionVector(element),...
%                     stress_projection_vector);
%                 mass_matrix = self.assembleStiffness(element, element.type.implementation.calcMassMatrix( element ),...
%                     mass_matrix);
%             end
%             nodal_stress_vector = mass_matrix\stress_projection_vector;
%             for node = self.model.node_dict
%                 node.setStressVonMises(nodal_stress_vector(node.number));
%             end
%         end        
    end
    
    methods( Static )
        function global_mass = assembleMassMatrix(node_container, local_mass, global_mass)
            %THIS IS THE SAME ASSEMBLER AS USED FOR THE STIFFNESS MATRIX IN
            %STANDARD SOOFEAM
            % dofs_per_node is dimension of stiffness matrix divided by number of nodes per element.
            % I have not figured out why, yet. Maybe number of entries
            % means numbers of dofs per node, so 2 for 2D and 3 for 3D
            % original line: nr_of_entries = local_stiffness.shape(1)/node_container.type.shape.getNumberOfNodes;
            
            dofs_per_element = length(local_mass);
            
            dofs_per_node = dofs_per_element/node_container.type.shape.getNumberOfNodes;
            
            for local_row_index = 1:length(local_mass)
                for local_col_index = 1:length(local_mass)
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
                    
                    global_mass(global_row_index,global_col_index) = global_mass(global_row_index,global_col_index) + local_mass(local_row_index,local_col_index);
                end
            end
        end

        function global_load = assembleLoad(node_container, local_load, global_load)
            
            % analogous to assembleStiffness(...) method above
               dofs_per_node = length(local_load)/node_container.type.shape.getNumberOfNodes;
            for local_row_index = 1:length(local_load)
                
                local_node_index_row = ceil(local_row_index/dofs_per_node);
                %row_node = node_container.node_list(local_node_index_row);
                %global_node_index_row = row_node.number;
                
                global_node_index_row = node_container.node_list(local_node_index_row).number;
                
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

%     methods ( Abstract )
%         run(self)
%     end
%         
%     methods ( Abstract, Access = protected )
%         updateDOF(self, solution_vector)
%         calcGlobalLoadVector(self, global_load)
%     end
end

