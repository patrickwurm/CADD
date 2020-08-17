function [pos, vel, xi, xi1, xi2, force, ekin, epot, temp, energies, virial] = Solve(pos, padpos, vel, xi, xi1, xi2, force, deltat, tau, mxatm, mxpatm, amass, rcut, xbox, ybox, tset, pbc, kb, nl, numneigh, aa, settings, example_name, force_tables)

sumv=0;
sumv2=0;

global inactive_atoms
global refpos
global exclude_atoms

thermostatted = settings.thermostatted;

atoms = [1:mxatm];
atoms(ismember(atoms,inactive_atoms))=[];
atoms(ismember(atoms,exclude_atoms))=[];

if thermostatted == 1
    % Algorithm: See "Thermostatting in Molecular Dynamics Simulations" pg. 7
    xi1= xi1 + xi2 * deltat/2;
    for ia=atoms
        vel(ia,:) = vel(ia,:) * exp( -xi1*deltat/2 ) + force(ia,:)/amass*deltat/2;
        pos(ia,:) = pos(ia,:) + vel(ia,:) * deltat;
    end
    %-------------------------------
    [force, en, energies, virial] = getforce(rcut, [pos; padpos], xbox, ybox, mxatm, mxpatm, pbc, nl, numneigh, settings, force_tables);
    %-------------------------------
    for ia=atoms
        vel(ia,:) = ( vel(ia,:) + force(ia,:)/amass*deltat/2 ) * exp( -xi1*deltat/2 );
        sumv2=sumv2+vel(ia,1)^2+vel(ia,2)^2;
    end
    temp = amass/(kb*(2*(mxatm-length(inactive_atoms)-length(exclude_atoms))-2))*sumv2;
    xi2 = 1/tau^2 * (temp/tset-1);
    xi1 = xi1 + xi2 * deltat/2;
    xi = xi + xi1*deltat;
    epot = en;
    ekin = 0.5*amass*sumv2;
    
else
    gamma_0 = settings.gamma_0;
    k = kb;
    
    if thermostatted == 2 || thermostatted == 3
        vel_hat = zeros(mxatm,2);
        a = zeros(mxatm,2);
        a_hat = zeros(mxatm,2);
        Fa = zeros(mxatm,2);
    end
    
    depth = settings.damping_band_depth;
    minx = min(refpos(:,1));
    miny = min(refpos(:,2));
    maxx = max(refpos(:,1));
    maxy = max(refpos(:,2));
    
    %Velocity Verlet
    for ia=atoms
        if strcmp(example_name,'Radial_pulse')
            if (thermostatted == 2 || thermostatted == 3) && (refpos(ia,2)<=miny+depth || refpos(ia,1) <= minx+depth || refpos(ia,1) >= maxx-depth || refpos(ia,2) >= maxy-depth)
                if thermostatted == 2
                    if refpos(ia,2)<=miny+depth
                        gamma = gamma_0/depth * (miny + depth - refpos(ia,2));
                    elseif refpos(ia,1) <= minx+depth
                        gamma = gamma_0/depth * (minx + depth - refpos(ia,1));
                    elseif refpos(ia,1) >= maxx-depth
                        gamma = gamma_0/depth * (refpos(ia,1) - (maxx - depth));
                    elseif refpos(ia,2) >= maxy-depth
                        gamma = gamma_0/depth * (refpos(ia,2) - (maxy - depth));
                    else
                        error('Something went wrong')
                    end
                elseif thermostatted == 3
                    gamma = gamma_0;
                end

                %square_distribution
                Fa(ia,:) = sqrt(6*gamma*amass*k*tset/deltat)*[(2*rand()-1) (2*rand()-1)];

                a(ia,:) = force(ia,:)/amass - gamma*vel(ia,:) + Fa(ia,:)/amass;
                pos(ia,:) = pos(ia,:) + vel(ia,:) * deltat + a(ia,:)*deltat^2/2;
            else
                pos(ia,:) = pos(ia,:) + vel(ia,:) * deltat + force(ia,:)/amass*deltat^2/2;
                vel(ia,:) = vel(ia,:) + force(ia,:)/amass*deltat/2; %half step velocity
            end
        elseif strcmp(example_name,'Example2') || strcmp(example_name,'Example2_rotated')
            if (thermostatted == 2 || thermostatted == 3)  && (refpos(ia,1) <= minx+depth)
                if thermostatted == 2
                    gamma = gamma_0/depth * (minx + depth - refpos(ia,1));
                elseif thermostatted == 3
                    gamma = gamma_0;
                end

                %square_distribution
                Fa(ia,:) = sqrt(6*gamma*amass*k*tset/deltat)*[(2*rand()-1) (2*rand()-1)];
                %+++++++++
                a(ia,:) = force(ia,:)/amass - gamma*vel(ia,:) + Fa(ia,:)/amass;
                pos(ia,:) = pos(ia,:) + vel(ia,:) * deltat + a(ia,:)*deltat^2/2;
            else
                pos(ia,:) = pos(ia,:) + vel(ia,:) * deltat + force(ia,:)/amass*deltat^2/2;
                vel(ia,:) = vel(ia,:) + force(ia,:)/amass*deltat/2; %half step velocity
            end
        elseif strcmp(example_name,'Tensile_Test') || strcmp(example_name,'Tensile_Test_rotated')
            if (thermostatted == 2 || thermostatted == 3)  && (refpos(ia,1) <= minx+depth || refpos(ia,1) >= maxx-depth)
                if thermostatted == 2
                    if refpos(ia,1) <= minx+depth
                        gamma = gamma_0/depth * (minx + depth - refpos(ia,1));
                    elseif refpos(ia,1) >= maxx-depth
                        gamma = gamma_0/depth * (refpos(ia,1) - (maxx - depth));
                    else
                        error('Something went wrong')
                    end
                elseif thermostatted == 3
                    gamma = gamma_0;
                end

                %square_distribution
                Fa(ia,:) = sqrt(6*gamma*amass*k*tset/deltat)*[(2*rand()-1) (2*rand()-1)];
                %+++++++++
                a(ia,:) = force(ia,:)/amass - gamma*vel(ia,:) + Fa(ia,:)/amass;
                pos(ia,:) = pos(ia,:) + vel(ia,:) * deltat + a(ia,:)*deltat^2/2;
            else
                pos(ia,:) = pos(ia,:) + vel(ia,:) * deltat + force(ia,:)/amass*deltat^2/2;
                vel(ia,:) = vel(ia,:) + force(ia,:)/amass*deltat/2; %half step velocity
            end
            
            
        else
            error('Example not implemented')
        end
    end
    %-------------------------------
    [force, en, energies, virial] = getforce(rcut, [pos; padpos], xbox, ybox, mxatm, mxpatm, pbc, nl, numneigh, settings, force_tables);
    %-------------------------------
    velcount = 0;
    for ia=atoms
        if strcmp(example_name,'Radial_pulse')
            if (thermostatted == 2 || thermostatted == 3)  && (refpos(ia,2)<=miny+depth || refpos(ia,1) <= minx+depth || refpos(ia,1) >= maxx-depth || refpos(ia,2) >= maxy-depth)
                %TEST STADIUM
                if thermostatted == 2
                    if refpos(ia,2)<=miny+depth
                        gamma = gamma_0/depth * (miny + depth - refpos(ia,2));
                    elseif refpos(ia,1) <= minx+depth
                        gamma = gamma_0/depth * (minx + depth - refpos(ia,1));
                    elseif refpos(ia,1) >= maxx-depth
                        gamma = gamma_0/depth * (refpos(ia,1) - (maxx - depth));
                    elseif refpos(ia,2) >= maxy-depth
                        gamma = gamma_0/depth * (refpos(ia,2) - (maxy - depth));
                    else
                        error('Something went wrong')
                    end
                elseif thermostatted == 3
                    gamma = gamma_0;
                end
                vel_hat(ia,:) = vel(ia,:) + 1/2*deltat*( a(ia,:) + force(ia,:)/amass - gamma*(vel(ia,:)+deltat*a(ia,:)) + Fa(ia,:)/amass  );
                a_hat(ia,:) = force(ia,:)/amass - gamma*(vel_hat(ia,:)) + Fa(ia,:)/amass;
                vel(ia,:) = vel(ia,:) + 1/2*deltat*(a(ia,:)+a_hat(ia,:));
            else
                vel(ia,:) = vel(ia,:) + force(ia,:)/amass*deltat/2;
                velcount = velcount + 1;
                sumv2=sumv2+vel(ia,1)^2+vel(ia,2)^2;
            end
        elseif strcmp(example_name,'Example2') || strcmp(example_name,'Example2_rotated')
            if (thermostatted == 2 || thermostatted == 3)  && refpos(ia,1) <= minx+depth
                if thermostatted == 2
                    gamma = gamma_0/depth * (minx + depth - refpos(ia,1));
                elseif thermostatted == 3
                    gamma = gamma_0;
                end
                vel_hat(ia,:) = vel(ia,:) + 1/2*deltat*( a(ia,:) + force(ia,:)/amass - gamma*(vel(ia,:)+deltat*a(ia,:)) + Fa(ia,:)/amass  );
                a_hat(ia,:) = force(ia,:)/amass - gamma*(vel_hat(ia,:)) + Fa(ia,:)/amass;
                vel(ia,:) = vel(ia,:) + 1/2*deltat*(a(ia,:)+a_hat(ia,:));
            else
                vel(ia,:) = vel(ia,:) + force(ia,:)/amass*deltat/2;
                velcount = velcount + 1;
                sumv2=sumv2+vel(ia,1)^2+vel(ia,2)^2;
            end
            
        elseif strcmp(example_name,'Tensile_Test') || strcmp(example_name,'Tensile_Test_rotated')
            if (thermostatted == 2 || thermostatted == 3)  && (refpos(ia,1) <= minx+depth || refpos(ia,1) >= maxx-depth)
                if thermostatted == 2
                    if refpos(ia,1) <= minx+depth
                        gamma = gamma_0/depth * (minx + depth - refpos(ia,1));
                    elseif refpos(ia,1) >= maxx-depth
                        gamma = gamma_0/depth * (refpos(ia,1) - (maxx - depth));
                    else
                        error('Something went wrong')
                    end
                elseif thermostatted == 3
                    gamma = gamma_0;
                end
                vel_hat(ia,:) = vel(ia,:) + 1/2*deltat*( a(ia,:) + force(ia,:)/amass - gamma*(vel(ia,:)+deltat*a(ia,:)) + Fa(ia,:)/amass  );
                a_hat(ia,:) = force(ia,:)/amass - gamma*(vel_hat(ia,:)) + Fa(ia,:)/amass;
                vel(ia,:) = vel(ia,:) + 1/2*deltat*(a(ia,:)+a_hat(ia,:));
            else
                vel(ia,:) = vel(ia,:) + force(ia,:)/amass*deltat/2;
                velcount = velcount + 1;
                sumv2=sumv2+vel(ia,1)^2+vel(ia,2)^2;
            end
            
        else
            error('Example not implemented')
        end
    end
    
    temp = amass/(kb*(2*velcount-2))*sumv2;
    
    epot = en;
    ekin = 0.5*amass*sumv2;
    
    xi=0;
    xi1=0;
    xi2=0;
end
end

