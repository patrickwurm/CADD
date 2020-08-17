function [force, en, energies, virial] = getforce(rcut, pos, xbox, ybox, mxatm, mxpatm, pbc, nl, numneigh, settings, force_tables)
en = 0;
force = zeros(mxatm,2);
energies = zeros(mxatm+mxpatm,1);

global indenter
global inactive_atoms
global virial_stress_atoms

%only implemented for EAM yet
virial = [];

if settings.atomic_potential == 1
    
    %%%%%%%%%%%%%%%%%%%%%
    %   EAM
    %%%%%%%%%%%%%%%%%%%%%
    
    rhoi = zeros(mxatm+mxpatm,1);
    Fi = zeros(mxatm+mxpatm,1);
    Fi1 = zeros(mxatm+mxpatm,1);
    
    if ~isempty(virial_stress_atoms)
    virial.sumvirial_xx = zeros(length(virial_stress_atoms),1);
    virial.sumvirial_yy = zeros(length(virial_stress_atoms),1);
    virial.sumvirial_xy = zeros(length(virial_stress_atoms),1);
    else
    virial.sumvirial_xx = 0;
    virial.sumvirial_yy = 0;
    virial.sumvirial_xy = 0;
    end
    virial_count = 0;
    
    convA_to_meters = 10^-10;
    conveV_to_J = 1.602176565*10^-19;
    
    %EAM Parameters
    %ZHOU 2001
    %re = 2.886166*convA_to_meters;
    %fe = 1.392302*conveV_to_J/convA_to_meters;
    %rhoe = 20.226537*conveV_to_J/convA_to_meters;
    %rhos = rhoe;
    %alpha = 6.942419;
    %beta = 3.702623;
    %A = 0.251519*conveV_to_J;
    %B = 0.313394*conveV_to_J;
    %kappa = 0.395132;
    %lambda = 0.790264;
    %Fn = [-2.806783; -0.276173; 0.893409; -1.637201]*conveV_to_J;
    %F = [-2.83; 0; 0.929508; -0.682320]*conveV_to_J;
    %eta = 0.779208;
    %Fe = -2.829437*conveV_to_J;

    %ZHOU 2004
    re = 2.863924*convA_to_meters;
    fe = 1.403115*conveV_to_J/convA_to_meters;
    rhoe = 20.418205*conveV_to_J/convA_to_meters;
    rhos = 23.195740*conveV_to_J/convA_to_meters;
    alpha = 6.613165;
    beta = 3.527021;
    %A = 0.134873*conveV_to_J; POSSIBLE MISTAKE IN ORIGINAL PAPER
    %https://atsimpotentials.readthedocs.io/en/latest/potentials/eam_tabulation.html#eam-example-2b
    %
    A = 0.314873*conveV_to_J; %
    B = 0.365551*conveV_to_J;
    kappa = 0.379846;
    lambda = 0.759692;
    Fn = [-2.807602; -0.301435; 1.258562; -1.247604]*conveV_to_J;
    F = [-2.83; 0; 0.622245; -2.488244]*conveV_to_J;
    eta = 0.785902;
    Fe = -2.824528*conveV_to_J;
    
    
    %First loop, calculate rhois and Fis
    %loop over all atoms (real + pad)
    for ia = 1:mxatm+mxpatm
        for ja=1:numneigh(ia)
            
            ka = nl(ia,ja);
            sx = pos(ia,1)-pos(ka,1);
            sy = pos(ia,2)-pos(ka,2);
            %PBC in y-direction (only valid in Example2)
            if pbc
                    if sy>ybox/2
                        sy=sy-ybox;
                    end
                    if sy<-ybox/2
                        sy=sy+ybox;
                    end
            end
            rrr = sqrt(sx^2 + sy^2);
            
            if rrr <= rcut
                %calculate rhois
                rhoi(ia) = rhoi(ia) +  fe*exp(-beta*(rrr/re-1))/ (1+(rrr/re-lambda)^20);
            end
        end
        
        %calculate Fis
        if rhoi(ia)<0.85*rhoe
            for k=0:3
                Fi(ia) = Fi(ia) + Fn(k+1)*(rhoi(ia)/(0.85*rhoe)-1)^k;
                Fi1(ia) = Fi1(ia) + Fn(k+1)*k*(rhoi(ia)/(0.85*rhoe)-1)^(k-1)*1/(0.85*rhoe);
            end
        elseif rhoi(ia) <1.15*rhoe
            for k=0:3
                Fi(ia) = Fi(ia) + F(k+1)*(rhoi(ia)/(rhoe)-1)^k;
                Fi1(ia) = Fi1(ia) + F(k+1)*k*(rhoi(ia)/(rhoe)-1)^(k-1)*1/rhoe;
            end
        else
            Fi(ia) = Fe*(1-log(rhoi(ia)/rhos)^eta)*(rhoi(ia)/rhos)^eta;
            Fi1(ia) = -Fe*eta/rhoi(ia)*(rhoi(ia)/rhos)^eta + Fe*(1-log(rhoi(ia)/rhos)^eta)*eta*(rhoi(ia)/rhos)^(eta-1)*1/rhos;
        end
    end
    
    sumphi=0;
    
    %force calculation
    %only loop over "real" atoms
    for ia = 1:mxatm
        %if isempty(find(inactive_atoms==ia))
        
        if ~isempty(find(virial_stress_atoms==ia,1))
            computevirial = 1;
            virial_count = virial_count + 1;
        else
            computevirial = 0;
        end
        
        vxx=0;
        vyy=0;
        vxy=0;
            
        for ja=1:numneigh(ia)
            
            ka = nl(ia,ja);
            sx = pos(ia,1)-pos(ka,1);
            sy = pos(ia,2)-pos(ka,2);
            
            if pbc
                    if sy>ybox/2
                        sy=sy-ybox;
                    end
                    if sy<-ybox/2
                        sy=sy+ybox;
                    end
            end
            
            rrr = sqrt(sx^2 + sy^2);
            if rrr <= rcut
                %calculate pi'
                pi1 = fe*exp(-beta*(rrr/re-1))*(-beta/re)*(1+(rrr/re-lambda)^20)^-1 ...
                    -fe*exp(-beta*(rrr/re-1))*(1+(rrr/re-lambda)^20)^-2*20*(rrr/re-lambda)^19*1/re;
                
                %calculate phi
                phi = A*exp(-alpha*(rrr/re-1))/(1+(rrr/re-kappa)^20)...
                    -B*exp(-beta*(rrr/re-1))/(1+(rrr/re-lambda)^20);
                
                %calculate phi'
                phi1 = A*exp(-alpha*(rrr/re-1))*(-alpha/re)*(1+(rrr/re-kappa)^20)^(-1) ...
                    - A*exp(-alpha*(rrr/re-1))*(1+(rrr/re-kappa)^20)^(-2)*20*(rrr/re-kappa)^19*1/re ...
                    - B*exp(-beta*(rrr/re-1))*(-beta/re)*(1+(rrr/re-lambda)^20)^(-1) ...
                    + B*exp(-beta*(rrr/re-1))*(1+(rrr/re-lambda)^20)^(-2)*20*(rrr/re-lambda)^19*1/re;
                
                factor = (phi1+Fi1(ia)*pi1+Fi1(ka)*pi1)/rrr;
                %factor*sx
                %factor*sy
                force(ia,1) = force(ia,1) - factor*sx;
                force(ia,2) = force(ia,2) - factor*sy;
                if computevirial
                    vxx=vxx-factor*sx*sx;
                    vyy=vyy-factor*sy*sy;
                    vxy=vxy-factor*sx*sy;
                    %virial.sumvirial_xx(virial_count) = virial.sumvirial_xx(virial_count) - factor*sx*sx;
                    %virial.sumvirial_yy(virial_count) = virial.sumvirial_yy(virial_count) - factor*sy*sy;
                    %virial.sumvirial_xy(virial_count) = virial.sumvirial_xy(virial_count) - factor*sx*sy;
                end
                sumphi = sumphi + phi;
                energies(ia) = energies(ia) + phi;
            end
            
            
        end

        if computevirial
            virial.sumvirial_xx(virial_count) = vxx;
            virial.sumvirial_yy(virial_count) = vyy;
            virial.sumvirial_xy(virial_count) = vxy;
        end
        en = en + Fi(ia);
        energies(ia) = energies(ia) + Fi(ia);
        %end
    end
    virial.sumvirial_xx = 1/2 * virial.sumvirial_xx;
    virial.sumvirial_yy = 1/2 * virial.sumvirial_yy;
    virial.sumvirial_xy = 1/2 * virial.sumvirial_xy;
    en = en + 1/2*sumphi;
    
elseif settings.atomic_potential == 2
    
    %%%%%%%%%%%%%%%%%%%%%
    %   LJ
    %%%%%%%%%%%%%%%%%%%%%
    
    %------------------------------------------------------------------
    %			real units	reduced units (ALUMINIUM)
    %------------------------------------------------------------------
    % epsilon = 4551*kb %epsilon in LJ Potential (see Halicioglu, T.; Pound, G. M. Phys. Status Solidi A 1975, 30.)
    % sigma = 2.62e-10 %sigma in LJ potential
    % ma = 26.981539*amu = 26.981539*1.66053904020e-27 || mass of aluminium
    
    epsilon = 4551 * 1.3806485279e-23;
    sigma = 2.62e-10;
    
    for ia = 1:mxatm
        if isempty(find(inactive_atoms==ia))
            for ja=1:numneigh(ia)
                ka = nl(ia,ja);
                sx = pos(ia,1)-pos(ka,1);
                sy = pos(ia,2)-pos(ka,2);
                %PBC in y-direction (only valid in Example2)
                if pbc
                    if sy>ybox/2
                        sy=sy-ybox;
                    end
                    if sy<-ybox/2
                        sy=sy+ybox;
                    end
                end
                
                rrr = sqrt(sx^2 + sy^2);
                
                if rrr <= rcut
                    %s6 = sigma^6;
                    %r6 = rrr^6;
                    
                    %factor = 24* s6 * epsilon * (r6 -2*s6)/(r6*rrr^7);
                    
                    factor = 24 * sigma^6 * epsilon * (rrr^6 - 2*sigma^6)/rrr^13;
                    
                    %force(ia,:) = force(ia,:) - factor/rrr*[sx,sy];
                    force(ia,1) = force(ia,1) - factor*sx/rrr;
                    force(ia,2) = force(ia,2) - factor*sy/rrr;
                    %en = en + 4*epsilon*( (sigma/rrr)^12 - (sigma/rrr)^6 ); %factor 1/2?
                    energies(ia) = energies(ia) + 4*( (sigma/rrr)^12 - (sigma/rrr)^6 );
                end
            end
        end
    end
    
    en = 0; %just to save time, not needed now
    
elseif settings.atomic_potential == 3
    %%%%%%%%%%%%%%%%%%%%%
    %   LJ (vectorized)
    %%%%%%%%%%%%%%%%%%%%%
    
    %------------------------------------------------------------------
    %			real units	reduced units (ALUMINIUM)
    %------------------------------------------------------------------
    % epsilon = 4551*kb %epsilon in LJ Potential (see Halicioglu, T.; Pound, G. M. Phys. Status Solidi A 1975, 30.)
    % sigma = 2.62e-10 %sigma in LJ potential
    % ma = 26.981539*amu = 26.981539*1.66053904020e-27 || mass of aluminium
    
    epsilon = 4551 * 1.3806485279e-23;
    sigma = 2.62e-10;
    
    %vectorized implementation of the LJ potential
    rrr_matrix = zeros([mxatm,size(nl,2)]);
    sxy_matrix = zeros([mxatm,size(nl,2),2]);
    
    for ia = 1:mxatm
        for ja=1:numneigh(ia)
            %ka = nl(ia,ja);
            sx = pos(ia,1)-pos(nl(ia,ja),1);
            sy = pos(ia,2)-pos(nl(ia,ja),2);
            %PBC in y-direction (only valid in Example2)
            if pbc
                    if sy>ybox/2
                        sy=sy-ybox;
                    end
                    if sy<-ybox/2
                        sy=sy+ybox;
                    end
            end

            sxy_matrix(ia,ja,1) = sx;
            sxy_matrix(ia,ja,2) = sy;
        end
    end
    % SET FORCES OF INACTIVE ATOMS TO ZERO (NOT NEEDED?)
    %sxy_matrix(inactive_atoms,:,:) = zeros(size(sxy_matrix(inactive_atoms,:,:)));
    
    rrr_matrix = sqrt(sxy_matrix(:,:,1).^2 + sxy_matrix(:,:,2).^2);
    rrr_matrix(rrr_matrix > rcut) = 0;
    
    factor_matrix = 24 * sigma^6 * epsilon * (rrr_matrix.^6 - 2*sigma^6)./rrr_matrix.^13;
    
    energy_matrix = 4*( (sigma./rrr_matrix).^12 - (sigma./rrr_matrix).^6 );
    
    factor_matrix(~isfinite(factor_matrix)) = 0;
    energy_matrix(~isfinite(energy_matrix)) = 0;
    
    force_matrix_1 = factor_matrix.*sxy_matrix(:,:,1)./rrr_matrix;
    force_matrix_2 = factor_matrix.*sxy_matrix(:,:,2)./rrr_matrix;
    force_matrix_1(~isfinite(force_matrix_1)) = 0;
    force_matrix_2(~isfinite(force_matrix_2)) = 0;
    
    force(:,1) = force(:,1) - sum(force_matrix_1,2);
    force(:,2) = force(:,2) - sum(force_matrix_2,2);
    
    energies = sum(energy_matrix,2);
    
    en = sum(energies);
    
    energies = vertcat(energies,zeros(mxpatm,1));
    
elseif settings.atomic_potential == 4
    
    %%%%%%%%%%%%%%%%%%%%%
    %   EAM (vectorized)
    %%%%%%%%%%%%%%%%%%%%%
    
    rhoi = zeros(mxatm+mxpatm,1);
    Fi = zeros(mxatm+mxpatm,1);
    Fi1 = zeros(mxatm+mxpatm,1);
    
    if ~isempty(virial_stress_atoms)
    virial.sumvirial_xx = zeros(length(virial_stress_atoms),1);
    virial.sumvirial_yy = zeros(length(virial_stress_atoms),1);
    virial.sumvirial_xy = zeros(length(virial_stress_atoms),1);
    else
    virial.sumvirial_xx = 0;
    virial.sumvirial_yy = 0;
    virial.sumvirial_xy = 0;
    end
    virial_count = 0;
    
    convA_to_meters = 10^-10;
    conveV_to_J = 1.602176565*10^-19;
    
    %EAM Parameters
    %ZHOU 2001
    %re = 2.886166*convA_to_meters;
    %fe = 1.392302*conveV_to_J/convA_to_meters;
    %rhoe = 20.226537*conveV_to_J/convA_to_meters;
    %rhos = rhoe;
    %alpha = 6.942419;
    %beta = 3.702623;
    %A = 0.251519*conveV_to_J;
    %B = 0.313394*conveV_to_J;
    %kappa = 0.395132;
    %lambda = 0.790264;
    %Fn = [-2.806783; -0.276173; 0.893409; -1.637201]*conveV_to_J;
    %F = [-2.83; 0; 0.929508; -0.682320]*conveV_to_J;
    %eta = 0.779208;
    %Fe = -2.829437*conveV_to_J;

    %ZHOU 2004
    re = 2.863924*convA_to_meters;
    fe = 1.403115*conveV_to_J/convA_to_meters;
    rhoe = 20.418205*conveV_to_J/convA_to_meters;
    rhos = 23.195740*conveV_to_J/convA_to_meters;
    alpha = 6.613165;
    beta = 3.527021;
    %A = 0.134873*conveV_to_J; POSSIBLE MISTAKE IN ORIGINAL PAPER
    %https://atsimpotentials.readthedocs.io/en/latest/potentials/eam_tabulation.html#eam-example-2b
    %
    A = 0.314873*conveV_to_J; %
    B = 0.365551*conveV_to_J;
    kappa = 0.379846;
    lambda = 0.759692;
    Fn = [-2.807602; -0.301435; 1.258562; -1.247604]*conveV_to_J;
    F = [-2.83; 0; 0.622245; -2.488244]*conveV_to_J;
    eta = 0.785902;
    Fe = -2.824528*conveV_to_J;
    
    size_nl = size(nl);
    cols = size_nl(2);
    
    %vectorized implementation of the LJ potential
    %rrr_matrix = zeros(size_nl);
    %sx_matrix = zeros(size_nl);
    %sy_matrix = zeros(size_nl);
    
    %VERY OLD VERSION
%     for ia = 1:mxatm+mxpatm
%         for ja=1:numneigh(ia)
%             %ka = nl(ia,ja);
%             sx = pos(ia,1)-pos(nl(ia,ja),1);
%             sy = pos(ia,2)-pos(nl(ia,ja),2);
%             %PBC in y-direction (only valid in Example2)
%             if pbc
%                 if sy>ybox/2
%                     sy=sy-ybox;
%                 end
%                 if sy<-ybox/2
%                     sy=sy+ybox;
%                 end
%             end
%             
%             sxy_matrix(ia,ja,1) = sx;
%             sxy_matrix(ia,ja,2) = sy;
%         end
%     end

%     for ia = 1:mxatm+mxpatm
%         nl_vec = nl(ia,:);
%         indices = nl_vec(nl_vec~=0);
%         tmp_vec = pos(indices,:);
%         vec = pos(ia,:).*ones(numneigh(ia),2) - tmp_vec;
%         if pbc
%             v2 = vec(:,2);
%             pbc_greater = v2 > ybox/2;
%             pbc_smaller = v2 < -ybox/2;
%             v2(pbc_greater) = v2(pbc_greater) -ybox;
%             v2(pbc_smaller) = v2(pbc_smaller) +ybox;
%             vec(:,2) = v2;
%         end
%         sxy_matrix(ia,1:numneigh(ia),:) = vec;
%     end

    pos_dummy = pos;
    pos_dummy(end+1,:) = 0;
    nl_dummy = nl;
    nl_dummy(nl_dummy==0) = length(pos_dummy(:,1));
    
    sx_matrix = cast(nl~=0,'double').*pos(:,1) - reshape(pos_dummy(nl_dummy,1),size(nl));
    sy_matrix = cast(nl~=0,'double').*pos(:,2) - reshape(pos_dummy(nl_dummy,2),size(nl));
    
    if pbc
        sy_matrix(sy_matrix > ybox/2) = sy_matrix(sy_matrix > ybox/2) - ybox;
        sy_matrix(sy_matrix < -ybox/2) = sy_matrix(sy_matrix < -ybox/2) + ybox;
    end
    
    
    rrr_matrix = sqrt(sx_matrix.^2 + sy_matrix.^2);
    % SET FORCES OF INACTIVE ATOMS TO ZERO (NOT NEEDED?)
    %sxy_matrix(inactive_atoms,:,:) = zeros(size(sxy_matrix(inactive_atoms,:,:)));
    
    rrr_matrix(rrr_matrix > rcut) = 0;
    
    sx_matrix_atoms = sx_matrix(1:mxatm,:);
    sy_matrix_atoms = sy_matrix(1:mxatm,:);
    
    rrr_matrix_atoms = rrr_matrix(1:mxatm,:);
    
    if force_tables.use_them && 0 %dont use it in this case, there is no increase in speed
        rhoi_matrix = zeros(size(rrr_matrix));
        for i = 1:cols
            rhoi_matrix(:,i) = qinterp1(force_tables.X, force_tables.rhoi, rrr_matrix(:,i));
        end
    else
        rhoi_matrix = fe*exp(-beta*(rrr_matrix/re-1))./ (1+(rrr_matrix/re-lambda).^20);
    end

    rhoi_matrix(rrr_matrix == 0) = 0;
    
    rhoi = sum(rhoi_matrix,2); % CHECK, IS IDENTICAL
    
    %NOT REALLY FEASIBLE TO VECTORIZE?
    for ia = 1:mxatm+mxpatm
        if rhoi(ia)<0.85*rhoe
            for k=0:3
                Fi(ia) = Fi(ia) + Fn(k+1)*(rhoi(ia)/(0.85*rhoe)-1)^k;
                Fi1(ia) = Fi1(ia) + Fn(k+1)*k*(rhoi(ia)/(0.85*rhoe)-1)^(k-1)*1/(0.85*rhoe);
            end
        elseif rhoi(ia) <1.15*rhoe
            for k=0:3
                Fi(ia) = Fi(ia) + F(k+1)*(rhoi(ia)/(rhoe)-1)^k;
                Fi1(ia) = Fi1(ia) + F(k+1)*k*(rhoi(ia)/(rhoe)-1)^(k-1)*1/rhoe;
            end
        else
            Fi(ia) = Fe*(1-log(rhoi(ia)/rhos)^eta)*(rhoi(ia)/rhos)^eta;
            Fi1(ia) = -Fe*eta/rhoi(ia)*(rhoi(ia)/rhos)^eta + Fe*(1-log(rhoi(ia)/rhos)^eta)*eta*(rhoi(ia)/rhos)^(eta-1)*1/rhos;
        end
    end

    if force_tables.use_them
        pi1_matrix = zeros(size(rrr_matrix_atoms));
        phi_matrix = zeros(size(rrr_matrix_atoms));
        phi1_matrix = zeros(size(rrr_matrix_atoms));
        for i = 1:cols
            pi1_matrix(:,i) = qinterp1(force_tables.X, force_tables.pi1, rrr_matrix_atoms(:,i));
            phi_matrix(:,i) = qinterp1(force_tables.X, force_tables.phi, rrr_matrix_atoms(:,i));
            phi1_matrix(:,i) = qinterp1(force_tables.X, force_tables.phi1, rrr_matrix_atoms(:,i));
        end
    else
        pi1_matrix = fe*exp(-beta*(rrr_matrix_atoms/re-1))*(-beta/re).*(1+(rrr_matrix_atoms/re-lambda).^20).^-1 ...
            -fe*exp(-beta*(rrr_matrix_atoms/re-1)).*(1+(rrr_matrix_atoms/re-lambda).^20).^-2*20.*(rrr_matrix_atoms/re-lambda).^19.*1/re;
   
        phi_matrix = A*exp(-alpha*(rrr_matrix_atoms/re-1))./(1+(rrr_matrix_atoms/re-kappa).^20)...
        -B*exp(-beta*(rrr_matrix_atoms/re-1))./(1+(rrr_matrix_atoms/re-lambda).^20);
    
        phi1_matrix = A*exp(-alpha*(rrr_matrix_atoms/re-1))*(-alpha/re).*(1+(rrr_matrix_atoms/re-kappa).^20).^(-1) ...
        - A*exp(-alpha*(rrr_matrix_atoms/re-1)).*(1+(rrr_matrix_atoms/re-kappa).^20).^(-2)*20.*(rrr_matrix_atoms/re-kappa).^19*1/re ...
        - B*exp(-beta*(rrr_matrix_atoms/re-1))*(-beta/re).*(1+(rrr_matrix_atoms/re-lambda).^20).^(-1) ...
        + B*exp(-beta*(rrr_matrix_atoms/re-1)).*(1+(rrr_matrix_atoms/re-lambda).^20).^(-2)*20.*(rrr_matrix_atoms/re-lambda).^19*1/re;
    end
      
    %TESTING
    if 0>1
    pi1_matrix_test = zeros(size(rrr_matrix_atoms));
    phi_matrix_test = zeros(size(rrr_matrix_atoms));
    phi1_matrix_test = zeros(size(rrr_matrix_atoms));
    for i = 1:cols
        pi1_matrix_test(:,i) = qinterp1(force_tables.X, force_tables.pi1, rrr_matrix_atoms(:,i));
        phi_matrix_test(:,i) = qinterp1(force_tables.X, force_tables.phi, rrr_matrix_atoms(:,i));
        phi1_matrix_test(:,i) = qinterp1(force_tables.X, force_tables.phi1, rrr_matrix_atoms(:,i));
    end
    
    error_pi1 = max(max(abs(pi1_matrix_test./pi1_matrix - 1)))*100 % in percent
    error_phi = max(max(abs(phi_matrix_test./phi_matrix - 1)))*100 % in percent
    error_phi1 = max(max(abs(phi1_matrix_test./phi1_matrix - 1)))*100 % in percent
    end
    
    %This is done in the interpolation function
    %pi1_matrix(rrr_matrix_atoms == 0) = 0;
    %phi_matrix(rrr_matrix_atoms == 0) = 0;
    %phi1_matrix(rrr_matrix_atoms == 0) = 0;
    
    %%%%%%%%%%%%%%%%
    
    
%     Fi1_neigh_matrix = zeros(size(pi1_matrix));
%     collength = length(pi1_matrix(1,:));
%     for ia = 1:mxatm
%         nl_vec = nl(ia,:);
%         indices = nl_vec(nl_vec~=0);
%         tmp_vec = Fi1(indices);
%         if length(tmp_vec) < collength
%             tmp_vec(collength) = 0; 
%         end
%         Fi1_neigh_matrix(ia,:) = tmp_vec';
%     end
    
    Fi1_dummy = Fi1;
    Fi1_dummy(end+1) = 0;
    nl_dummy = nl(1:mxatm,:);
    nl_dummy(nl_dummy==0) = length(Fi1_dummy);
    
    Fi1_neigh_matrix = reshape(Fi1_dummy(nl_dummy),size(nl_dummy));

    factor_matrix = (phi1_matrix + Fi1(1:mxatm).*pi1_matrix + Fi1_neigh_matrix.*pi1_matrix)./rrr_matrix_atoms;
    
    factor_matrix(~isfinite(factor_matrix)) = 0;
    
    %force_matrix_1 = factor_matrix.*sx_matrix_atoms;
    %force_matrix_2 = factor_matrix.*sy_matrix_atoms;

    force(:,1) = - sum(factor_matrix.*sx_matrix_atoms,2);
    force(:,2) = - sum(factor_matrix.*sy_matrix_atoms,2);
    
    sum_phi_vec = sum(phi_matrix,2);
    sum_phi_vec(mxatm+mxpatm) = 0;
    energies = energies + sum_phi_vec +Fi;
    
    en = 1/2*(sum(sum(phi_matrix)));
    
    virial.sumvirial_xx = -1/2*sum(factor_matrix.*sx_matrix_atoms.*sx_matrix_atoms,2);
    virial.sumvirial_yy = -1/2*sum(factor_matrix.*sx_matrix_atoms.*sy_matrix_atoms,2);
    virial.sumvirial_xy = -1/2*sum(factor_matrix.*sx_matrix_atoms.*sy_matrix_atoms,2);
end

en = en*1/2; %factor 1/2?

end
