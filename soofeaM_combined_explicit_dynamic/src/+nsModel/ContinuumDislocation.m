classdef ContinuumDislocation < nsModel.NumberedObject
    %CONTINUUMDISLOCATION
    properties%( SetAccess = protected )
        coordinates
        direction %i guess we need this
        residing_element_number
        intensity % this variable (goes from 0 to 1 and) is used to gradually include the effects of the dislocation on the system
        %dislocation to the problem (this is needed to not have a sudden load on the dynamical problem (it will explode))
        birth_time % at which timestep was the dislocation born
        grow_up_time % how many timesteps does it take for the dislocation to reach its full potential
        drag_coefficient % drag coefficient B 
        delta_t_D %update timestep for discrete dislocations
        stress_DD 
        FE_stress_array
    end
    
    methods
        function self = ContinuumDislocation( number, coordinates, direction, element_dict, birth_time, grow_up_time, drag_coeff, delta_t_D, irunDD )
            self@nsModel.NumberedObject(number);
            self.coordinates = coordinates;
            self.direction = direction;
            self.residing_element_number = self.findElement(element_dict);
            self.intensity = 0;
            self.birth_time = birth_time;
            self.grow_up_time = grow_up_time;
            self.drag_coefficient = drag_coeff;
            self.delta_t_D = delta_t_D;
            self.stress_DD = [0;0;0];
            self.FE_stress_array = zeros(irunDD,1);
        end
        
        %This is pretty much the same function that is used to find the
        %elements in which the pad atoms reside
        function element_number = findElement(self, element_dict)
            associated = false;
            for element = element_dict
                if self.triangleTest(element.node_list, self.coordinates)
                    element_number = element.number;
                    associated = true;
                end
            end
            if associated == false
                error(['Continuum dislocation ', num2str(self.number),' is not associated with any element.'])
            end
        end
        
        function  in_triangle  = triangleTest(self, nodes, position_PA )
            %Test if a point with coordinates "position_PA" lies in the triangle
            %formed by the 3 nodes "nodes".
            
            a = [nodes(1).material_coordinates, 0];
            b = [nodes(2).material_coordinates, 0];
            c = [nodes(3).material_coordinates, 0];
            
            global tol
            criteria = tol*1e-35; %the tolerance here must be very tight to ensure proper association of the pad atoms
            
            criteria = tol*1e-34;
            
            p = [position_PA, 0];
            
            % activate this for a visualization of the pad atoms and the elements to
            % be tested
            if 0>1
                figure(71)
                line([a(1),b(1)], [a(2), b(2)], 'Color', 'r');
                hold on;
                line([b(1),c(1)], [b(2), c(2)], 'Color', 'r');
                line([a(1),c(1)], [a(2), c(2)], 'Color', 'r');
                plot(p(1),p(2),'k^');
                drawnow
            end
            
            if sameSide(p, a, b, c) && sameSide(p, b, a, c) && sameSide (p, c, a, b)
                in_triangle = true;
            else
                in_triangle = false;
            end
            
            function on_same_side = sameSide(p1, p2, a, b)
                cp1 = cross(b-a, p1-a);
                cp2 = cross(b-a, p2-a);
                if dot(cp1,cp2) >= -criteria
                    on_same_side = true;
                else
                    on_same_side = false;
                end
            end
        end
        
        function upDatePosition(self, element, a)
            if self.direction == 1 %CHECK IF THIS IS RIGHT
                b = a;
            else
                b = -a;
            end
            %OLD, NO AVERAGING OF FE STRESS:
            %shear_stress = element.int_points(1).stress(3) +
            %self.stress_DD(3);
            
            %NEW:
            shear_stress = mean(self.FE_stress_array) + self.stress_DD(3);
            
            disp(['Dislocation ',num2str(self.number),' in element ',num2str(self.residing_element_number), ' : S12_FE: ',num2str(mean(self.FE_stress_array)/10e6),' S12_DD: ',num2str(self.stress_DD(3)/10e6)])
            %self.coordinates(2) = self.coordinates(2) + self.drag_coefficient * shear_stress * abs(b) * self.delta_t_D;
            self.coordinates(2) = self.coordinates(2) + 1/self.drag_coefficient * shear_stress * b * self.delta_t_D;
        end
        
        function displacement = getDisplacementAtAtom(self, atom_coords, a, material)
            if self.direction == 1 %CHECK IF THIS IS RIGHT
                b = a;
            else
                b = -a;
            end
            
            % we gotta get Poissons number from the FE model
            % (somehow)
            nu = material.nu; %Poissons number
            
            dx2 = atom_coords(1) - self.coordinates(1);
            dx1 = atom_coords(2) - self.coordinates(2);
            
            if abs(dx1) < a/100 && abs(dx2) < a/100;
                displacement = 0;
            else
                % Needleman: wrong??
                displacement_y = b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) );
                displacement_x = b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) );
                displacement = [displacement_x, displacement_y];
                displacement = displacement*self.intensity; 
            end
        end
        
        function incremental_displacement = getIncrementalDisplacementAtAtom(self, atom_coords, a, material)
            if self.direction == 1 %CHECK IF THIS IS RIGHT
                b = a;
            else
                b = -a;
            end
            
            % we gotta get Poissons number from the FE model
            % (somehow)
            nu = material.nu; %Poissons number
            
            dx2 = atom_coords(1) - self.coordinates(1);
            dx1 = atom_coords(2) - self.coordinates(2);
            
            if abs(dx1) < a/100 && abs(dx2) < a/100;
                incremental_displacement = 0;
            else
                % Needleman: wrong??
                displacement_y = b/(2*pi*(1-nu)) * ( 1/2*(dx1.*dx2)/(dx1.^2+dx2.^2) - (1-nu)*atan(dx1./dx2) );
                displacement_x = b/(2*pi*(1-nu)) * ( 1/2*(dx2.^2)/(dx1.^2+dx2.^2) - 1/4*(1-2*nu)*log((dx1.^2+dx2.^2)/b^2) );
                displacement = [displacement_x, displacement_y];
                incremental_displacement = displacement/self.grow_up_time;
            end
        end
        
        function intensity = computeIntensity(self, timestep)
            if timestep - self.birth_time > self.grow_up_time
                self.intensity = 1;
                intensity = 1;
            else
                self.intensity = (timestep - self.birth_time)/self.grow_up_time;
                intensity = self.intensity;
            end
            
%             if self.number == 1 || self.number == 2 
%                 if timestep > 1500
%                     intensity = 0;
%                 end
%             end
        end
        
        function setStressDD(self, stress_DD)
            self.stress_DD = stress_DD;
        end
        
                
        function setResidingElement(self, element_number)
            self.residing_element_number = element_number;
        end
        
        function addFEStressToArray(self, stress_FE)
            self.FE_stress_array(1:end-1) = self.FE_stress_array(2:end);
            self.FE_stress_array(end) = stress_FE;
        end
        
        
    end
    
   
    
    
    
    
end