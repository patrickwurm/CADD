classdef Linear_Dynamic_Analysis < nsAnalyzer.nsAnalysis.Explicit_Dynamic_Analysis
    methods
        function self = Linear_Dynamic_Analysis(model, output_handler, example_name, include_external_BCs)
            if nargin > 3
            else
                include_external_BCs = 1;
            end
            self@nsAnalyzer.nsAnalysis.Explicit_Dynamic_Analysis(model, output_handler, example_name, include_external_BCs);
        end
        
        function paddisp = run(self, write_to_nodes, boundary_cond_matrix_interface, fe_step, ref_positions, positions, mxatm, mxpatm, updateDD, istep, settings)
            timer_1 = tic;
            i_percent = 1;
            
            % Read (FE-)BC from input file
            self.boundary_cond_matrix = self.model.bc_handler.incorporateBC(self.model, 1, istep, settings);
            
            if self.include_external_BCs == 0
                self.boundary_cond_matrix(:,2:4) = 0;
            end
            
            % APPEND MD BC to boundary_cond_matrix
            %self.boundary_cond_matrix = vertcat(self.boundary_cond_matrix, boundary_cond_matrix_interface);
            self.boundary_cond_matrix = vertcat(boundary_cond_matrix_interface, self.boundary_cond_matrix);
            
            %TEST ADD NEW CONT DISLOCATION (DONT FORGET TO DELETE ME LATER)
%             if isempty(self.model.dislocation_dict)
%                 self.model.addContinuumDislocation(nsModel.ContinuumDislocation(length(self.model.dislocation_dict)+1, [15*self.model.a, -2*self.model.a], 1, self.model.element_dict));
%             end
            
            % Incorporate displacements from discrete dislocations (real and ghost)
            if ~isempty(self.model.dislocation_dict)% check if there are any continuum dislocations
                % first call, write nodal coordinates into matrix
                % This wont change during the run
                if isempty(self.DD_nodal_coordinates)
                    for i = 1:length(self.boundary_cond_matrix(:,5))
                        self.DD_nodal_coordinates(i,:) = self.model.node_dict(self.boundary_cond_matrix(i,5)).material_coordinates;
                    end
                end
                %%%%%%%%%%%
                %CARE
                %%%%%%%%%%%
                self.boundary_cond_matrix = self.includeDDDisplacements;
            end
            
            self.solveFESystem;
           
            if updateDD
                write_to_nodes = 1;
            end
            %write displacements to nodes only when needed for the output
            if write_to_nodes
                %compute u_dd (this is enough since we dont really need it every timestep?)
                if ~isempty(self.model.dislocation_dict)
                    for node = self.model.node_dict
                        for coord_index = 1:self.model.dimension
                            global_index = (node.number-1) * self.model.dimension + coord_index;
                            self.u_dd(global_index) = self.getDDDisplacementatNode(node, coord_index);
                        end
                    end
                end
                self.updateDOF(self.u + self.u_dd, self.u)
            end
            
            %compute stresses at needed elements
            if ~isempty(self.model.dislocation_dict)
                self.updateDOF(self.u + self.u_dd, self.u) %we need u (not u_dd) to be up_to_date in the nodes
                for i = self.model.dislocation_dict
                    element = self.model.element_dict(i.residing_element_number);
                    element.type.implementation.calcStrainStressInIP(element);
                    %and add the stress to the array of stresses in the
                    %dislocation (and average later)
                    i.addFEStressToArray(element.int_points(1).stress(3));
                end
            end
                
            if updateDD
                %compute stresses in all Gauss points of elements which
                %have dislocations in them
                if ~isempty(self.model.dislocation_dict)
                    %%%%%%%%%%%
                    %CARE
                    %%%%%%%%%%%
                    self.computeDDStressAtDislocationSites;
                end
                
                %Update dislocation positions and
                %find the element, theyre in now
                for dislocation = self.model.dislocation_dict
                    dislocation.upDatePosition(self.model.element_dict(dislocation.residing_element_number),self.model.a);
                    %dislocation.findElement(self.model.element_dict); %not
                    %needed? always take stress at "material position" and
                    %dont update the position?
                end
            end
            
            %Get Pad atom displacements
            paddisp = self.getPadDisplacements;
            
            %Update Continuum Dislocation Positions
            
            %self.calcStrainStressForPostprocessing()
        end
        
    end
    
    methods
        
        function updateDOF(self, solution_vector, solution_vector_FE)
            for node = self.model.node_dict
                for coord_index = 1:self.model.dimension
                    global_index = (node.number-1) * self.model.dimension + coord_index;
                    node.dof.setDisplacement(coord_index, solution_vector(global_index));
                    node.dof.setDisplacementFE(coord_index, solution_vector_FE(global_index));
                end
                
            end
        end
        
    end
    
    methods (Access = protected)    
        
        function incorporateAtomicDisplacements(self, interface_atom_displacements)
            for node = self.model.node_dict
                for i = 1:length(interface_atom_displacements(:,1))
                    if node.interface_atom_number == interface_atom_displacements(i,1)
                        node.setBCDOF('x',interface_atom_displacements(i,2));
                        node.setBCDOF('y',interface_atom_displacements(i,3));
                    end
                end
            end
        end
        
        function incorporateAtomicForces(self, interface_atom_forces)
            for node = self.model.node_dict
                for i = 1:length(interface_atom_displacements(:,1))
                    if node.interface_atom_number == interface_atom_displacements(i,1)
                        node.setNodeForce(interface_atom_forces(i,:));
                    end
                end
            end
        end
        
%         function padpos = updatePadAtoms(self, mxatm, mxpatm)
%             padpos = zeros(mxpatm,2);
%             
%             for element = self.model.element_dict
%                 for pad_atom = element.pad_atoms
%                     displacements = element.getDisplacementArray();
%                     H = element.type.shape.getArray( pad_atom.natural_coordinates );
%                     pad_atom.setDisplacements( H*displacements );
%                     padpos(pad_atom.number-mxatm,:) = pad_atom.coordinates + pad_atom.displacements;
%                 end
%             end
%         end

        
        function global_load = calcGlobalLoadVector(self, global_load)
            %disp('\--> Assemble global load')
            
            % Internal Forces:
            for element = self.model.element_dict
                F_int = element.type.implementation.calcLoad( element );
                global_load = self.assembleLoad(element, F_int, global_load);
            end
            
            % NOT NEEDED (no external forces in CADD)
            %global_load = calcExternalLoadVector(self, global_load);
        end
        
        function updated_bc_matrix = includeDDDisplacements(self)
            
            nu = self.model.material_dict{1,1}.nu; %Poissons number
            
            bc_matrix = self.boundary_cond_matrix;
            for i = 1:length(bc_matrix(:,1))
                % 1) LOOP OVER ALL CONTINUUM DISLOCATIONS
                for j = self.model.dislocation_dict
                    if j.direction == 1 %CHECK IF THIS IS RIGHT
                        b = self.model.a;
                    else
                        b = -self.model.a;
                    end
                    
                    dx2 = self.DD_nodal_coordinates(i,1) - j.coordinates(1);
                    dx1 = self.DD_nodal_coordinates(i,2) - j.coordinates(2);
                    
                    if abs(dx1) < self.model.a/100 && abs(dx2) < self.model.a/100;
                        %skip
                    else
                        if bc_matrix(i,6) == 1 %x-direction
                            bc_matrix(i,2) = bc_matrix(i,2) - b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) )*j.computeIntensity(self.model.time);
                        else %y-direction
                            bc_matrix(i,2) = bc_matrix(i,2) - b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) )*j.computeIntensity(self.model.time);
                        end
                    end
                end
                %2) LOOP OVER ALL GHOST DISLOCATIONS
                for j = self.model.ghostdislocation_dict
                    if j.direction == 1 %CHECK IF THIS IS RIGHT
                        b = self.model.a;
                    else
                        b = -self.model.a;
                    end
                    
                    dx2 = self.DD_nodal_coordinates(i,1) - j.coordinates(1);
                    dx1 = self.DD_nodal_coordinates(i,2) - j.coordinates(2);
                    
                    if abs(dx1) < self.model.a/100 && abs(dx2) < self.model.a/100;
                        %skip
                    else
                        if bc_matrix(i,6) == 1 %x-direction
                            bc_matrix(i,2) = bc_matrix(i,2) - b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) )*j.computeIntensity(self.model.time);
                        else %y-direction
                            bc_matrix(i,2) = bc_matrix(i,2) - b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) )*j.computeIntensity(self.model.time);
                        end
                    end
                    
                end
            end
            updated_bc_matrix = bc_matrix;
        end
        
        function displacement = getDDDisplacementatNode(self, node, coord_index)
            
            nu = self.model.material_dict{1,1}.nu; %Poissons number
            displacement = 0;
            % 1) LOOP OVER ALL CONTINUUM DISLOCATIONS
            for j = self.model.dislocation_dict
                if j.direction == 1 %CHECK IF THIS IS RIGHT
                    b = self.model.a;
                else
                    b = -self.model.a;
                end
                
                dx2 = node.material_coordinates(1) - j.coordinates(1);
                dx1 = node.material_coordinates(2) - j.coordinates(2);
                
                if abs(dx1) < self.model.a/100 && abs(dx2) < self.model.a/100;
                    displacement = 0;
                else
                    if coord_index == 1 %x-direction
                        displacement = displacement + b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) )*j.computeIntensity(self.model.time);;
                    else %y-direction
                        displacement = displacement + b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) )*j.computeIntensity(self.model.time);;
                    end
                end
            end
            % 2) LOOP OVER ALL GHOST DISLOCATIONS
            for j = self.model.ghostdislocation_dict
                if j.direction == 1 %CHECK IF THIS IS RIGHT
                    b = self.model.a;
                else
                    b = -self.model.a;
                end
                
                dx2 = node.material_coordinates(1) - j.coordinates(1);
                dx1 = node.material_coordinates(2) - j.coordinates(2);
                
                if abs(dx1) < self.model.a/100 && abs(dx2) < self.model.a/100;
                    displacement = 0;
                else
                    if coord_index == 1 %x-direction
                        displacement = displacement + b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) )*j.computeIntensity(self.model.time);
                    else %y-direction
                        displacement = displacement + b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) )*j.computeIntensity(self.model.time);
                    end
                end
            end
            
        end
        
        function computeDDStressAtDislocationSites(self)
            
            % we gotta get Poissons number and shear modulus from the FE model
            % (somehow)
            nu = self.model.material_dict{1,1}.nu; %Poissons number
            mu = self.model.material_dict{1,1}.mu; %Shear modulus
            for k = self.model.dislocation_dict
                s11 = 0;
                s22 = 0;
                s12 = 0;
                % 1) LOOP OVER ALL CONTINUUM DISLOCATIONS (but not over itself)
                for j = self.model.dislocation_dict
                    if j ~= k
                        if j.direction == 1 %CHECK IF THIS IS RIGHT
                            b = self.model.a;
                        else
                            b = -self.model.a;
                        end
                        
                        dx2 = k.coordinates(1) - j.coordinates(1);
                        dx1 = k.coordinates(2) - j.coordinates(2);
                        
                        if abs(dx1) < self.model.a/100 && abs(dx2) < self.model.a/100;
                            %do nothing
                        else
                            s22 = s22 - mu*b/(2*pi*(1-nu))*(dx2*(3*dx1^2 + dx2^2))/(dx1^2 + dx2^2)^2*j.computeIntensity(self.model.time);
                            s11 = s11 + mu*b/(2*pi*(1-nu))*(dx2*(dx1^2 - dx2^2))/(dx1^2 + dx2^2)^2*j.computeIntensity(self.model.time);
                            s12 = s12 + mu*b/(2*pi*(1-nu))*(dx1*(dx1^2 - dx2^2))/(dx1^2 + dx2^2)^2*j.computeIntensity(self.model.time);
                        end
                    end
                end
                % 2) LOOP OVER ALL GHOST DISLOCATIONS
                for j = self.model.ghostdislocation_dict
                    if j.direction == 1 %CHECK IF THIS IS RIGHT
                        b = self.model.a;
                    else
                        b = -self.model.a;
                    end
                    
                    % we gotta get Poissons number and shear modulus from the FE model
                    % (somehow)
                    nu = self.model.material_dict{1,1}.nu; %Poissons number
                    mu = self.model.material_dict{1,1}.mu; %Shear modulus
                    
                    dx2 = k.coordinates(1) - j.coordinates(1);
                    dx1 = k.coordinates(2) - j.coordinates(2);
                    
                    if abs(dx1) < self.model.a/100 && abs(dx2) < self.model.a/100;
                        %do nothing
                    else
                        s22 = s22 - mu*b/(2*pi*(1-nu))*(dx2*(3*dx1^2 + dx2^2))/(dx1^2 + dx2^2)^2*j.computeIntensity(self.model.time);
                        s11 = s11 + mu*b/(2*pi*(1-nu))*(dx2*(dx1^2 - dx2^2))/(dx1^2 + dx2^2)^2*j.computeIntensity(self.model.time);
                        s12 = s12 + mu*b/(2*pi*(1-nu))*(dx1*(dx1^2 - dx2^2))/(dx1^2 + dx2^2)^2*j.computeIntensity(self.model.time);
                    end
                end
                k.setStressDD([s11;s22;s12]);
            end
        end

    end
end

