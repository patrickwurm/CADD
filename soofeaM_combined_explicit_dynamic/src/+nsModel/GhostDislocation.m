classdef GhostDislocation < nsModel.NumberedObject
    %CONTINUUMDISLOCATION
    properties( SetAccess = protected )
        coordinates
        direction %i guess we need thisd
        intensity % this variable (goes from 0 to 1 and) is used to gradually include the effects of the dislocation on the system        %the following 2 variables are useful to gradually introduce the
        %dislocation to the problem (this is needed to not have a sudden load on the dynamical problem (it will explode))
        birth_time % at which timestep was the dislocation born
        grow_up_time % how many timesteps does it take for the dislocation to reach its full potential
    end
    
    methods
        function self = GhostDislocation( number, coordinates, direction, birth_time, grow_up_time )
            self@nsModel.NumberedObject(number);
            self.coordinates = coordinates;
            self.direction = direction;
            self.intensity = 0;
            self.birth_time = birth_time;
            self.grow_up_time = grow_up_time;
        end
        
        
        function displacement = getGhostDisplacementAtAtom(self, atom_coords, a, material)
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
        
        function incremental_displacement = getGhostIncrementalDisplacementAtAtom(self, atom_coords, a, material)
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
        end
        
    end
    
end