function [pbc_atoms, band_atoms, band_atoms_pos, band_atom_ass_atoms, band_atom_ass_atoms_position, inactive_atoms, exclude_atoms, lower_edge_atoms, elements, virial_stress_atoms] = getSpecialAtoms(example_name, pos, refpos, refpadpos, xbox, ybox, bandatom_position, a, axmin, axmax, aymin, aymax, rcut, inactive_atoms, exclude_atoms, fe_model, interface_atom_list, test_band_atom_association,dx , dy)
% This function outputs all "special" atoms and is very problem specific

pbc_atoms = [];
global tol

if strcmp(example_name, 'Example2')
    band_atoms = [];
    band_atoms_mean = [];
    lower_edge_atoms = [];
    pbc_atoms = []; %list of ''excluded'' atoms (cause of PBC), care must be taken with these
    for i=1:length(pos) % ONLY real atoms, no pad atoms, this is different in neighborlist
        if abs(pos(i,2)-ybox)<=10^-15 %careful, comparison with ybox is only valid if ymin = 0!
            for j = 1:length(pos)
                if abs(pos(j,2)-0)<=10^-15 && abs(pos(i,1)-pos(j,1))<=10^-15
                    pbc_atoms = [pbc_atoms, [i;j]];
                end
            end
        end
        if abs(pos(i,1)-bandatom_position*dx)<=10^-15
            band_atoms = [band_atoms, i];
        end
    end
    band_atoms(end) = []; % The upmost band atom is equal to the lowest one!
    band_atoms_mean = zeros(length(band_atoms),2);
    band_atoms_pos = zeros(size(band_atoms_mean));
    for i = 1:length(band_atoms)
        band_atoms_pos(i,:) = refpos(band_atoms(i),:);
        band_atoms_mean(i,:) = refpos(band_atoms(i),:);
    end
    
    xbandl = round(bandatom_position/(dx/a))*dx;
    
    %Associate a node to every band atom (needed as a reference for the new algorithm)
    disp('Associating Band Atoms with Interface Atoms ...')
    band_atom_ass_atoms = associateIntAtomswBandAtoms(example_name, band_atoms, pos, fe_model,  xbandl, 0, 0, 0, interface_atom_list);
    band_atom_ass_atoms_position = zeros(length(band_atom_ass_atoms),2);
    
elseif strcmp(example_name, 'Dislocation')
    bandatom_position = settings.bandatom_position;
    band_atoms = [];
    band_atoms_mean = [];
    
    %this all just works if xmin, ymin = 0
    yband = bandatom_position*dy;
    xbandl = round(bandatom_position/(dx/a))*dx;
    xbandr = xbox - xbandl;
    
    for i=1:length(pos) % ONLY real atoms, no pad atoms, this is different in neighborlist
        if abs(pos(i,2)-yband)<=10^-15 && pos(i,1) >= xbandl - 10^-15 && pos(i,1) <= xbandr + 10^-15
            band_atoms(end+1) = i;
        end
        
        if abs(pos(i,1)-xbandl)<=10^-15 && pos(i,2) >= yband - 10^-15
            band_atoms(end+1) = i;
        end
        
        if abs(pos(i,1)-xbandr)<=10^-15 && pos(i,2) >= yband - 10^-15
            band_atoms(end+1) = i;
        end
    end
    band_atoms = unique(band_atoms);
    band_atoms_mean = zeros(length(band_atoms),2);
    band_atoms_pos = zeros(size(band_atoms_mean));
    
    
    %Associate a node to every band atom (needed as a reference for the new algorithm)
    disp('Associating Band Atoms with Interface Atoms ...')
    band_atom_ass_atoms = associateIntAtomswBandAtoms(example_name, band_atoms, pos, fe_model, xbandl, xbandr, yband, interface_atom_list);
    band_atom_ass_atoms_position = zeros(length(band_atom_ass_atoms),2);
    
elseif strcmp(example_name, 'Radial_Pulse')
    band_atoms = [];
    band_atoms_mean = [];
    
    %this all just works if xmin, ymin = 0
    ybandb = bandatom_position*dy;
    ybandt = ybox-ybandb;
    xbandl = round(bandatom_position/(dx/a))*dx;
    xbandr = xbox - xbandl;
    
    for i=1:length(pos) % ONLY real atoms, no pad atoms, this is different in neighborlist
        if abs(pos(i,2)-ybandb)<=10^-15 && pos(i,1) >= xbandl - 10^-15 && pos(i,1) <= xbandr + 10^-15
            band_atoms(end+1) = i;
        end
        
        if abs(pos(i,2)-ybandt)<=10^-15 && pos(i,1) >= xbandl - 10^-15 && pos(i,1) <= xbandr + 10^-15
            band_atoms(end+1) = i;
        end
        
        if abs(pos(i,1)-xbandl)<=10^-15 && pos(i,2) >= ybandb - 10^-15 && pos(i,2) <= ybandt + 10^-15
            band_atoms(end+1) = i;
        end
        
        if abs(pos(i,1)-xbandr)<=10^-15 && pos(i,2) >= ybandb - 10^-15 && pos(i,2) <= ybandt + 10^-15
            band_atoms(end+1) = i;
        end
    end
    band_atoms = unique(band_atoms);
    band_atoms_mean = zeros(length(band_atoms),2);
    band_atoms_pos = zeros(size(band_atoms_mean));
   
    %Associate a node to every band atom (needed as a reference for the new algorithm)
    disp('Associating Band Atoms with Interface Atoms ...')
    band_atom_ass_atoms = associateIntAtomswBandAtoms(example_name, band_atoms, pos, fe_model, xbandl, xbandr, ybandb, ybandt, interface_atom_list);
    band_atom_ass_atoms_position = zeros(length(band_atom_ass_atoms),2);
   
elseif strcmp(example_name, 'Example2_rotated')
    band_atoms = [];
    band_atoms_mean = [];
    lower_edge_atoms = [];
    pbc_atoms = []; %list of ''excluded'' atoms (cause of PBC), care must be taken with these
    for i=1:length(pos) % ONLY real atoms, no pad atoms, this is different in neighborlist
        if abs(pos(i,2)-ybox)<=10^-15 %careful, comparison with ybox is only valid if ymin = 0!
            for j = 1:length(pos)
                if abs(pos(j,2)-0)<=10^-15 && abs(pos(i,1)-pos(j,1))<=10^-15
                    pbc_atoms = [pbc_atoms, [i;j]];
                end
            end
        end
        if abs(pos(i,1)-bandatom_position*dx)<=10^-15
            band_atoms = [band_atoms, i];
        end
    end
    band_atoms(end) = []; % The upmost band atom is equal to the lowest one!
    band_atoms_mean = zeros(length(band_atoms),2);
    band_atoms_pos = zeros(size(band_atoms_mean));
    for i = 1:length(band_atoms)
        band_atoms_pos(i,:) = refpos(band_atoms(i),:);
        band_atoms_mean(i,:) = refpos(band_atoms(i),:);
    end
    
    xbandl = round(bandatom_position)*dx;
    
    %Associate a node to every band atom (needed as a reference for the new algorithm)
    disp('Associating Band Atoms with Interface Atoms ...')
    band_atom_ass_atoms = associateIntAtomswBandAtoms(example_name, band_atoms, pos, fe_model,  xbandl, 0, 0, 0, interface_atom_list);
    band_atom_ass_atoms_position = zeros(length(band_atom_ass_atoms),2);

elseif strcmp(example_name, 'Tensile_Test')
    band_atoms = [];
    band_atoms_mean = [];
    lower_edge_atoms = [];
    pbc_atoms = [];
    
    xbandl = round(bandatom_position/(dx/a))*dx;
    xbandr = xbox - xbandl;
    
    pbc_atoms = []; %list of ''excluded'' atoms (cause of PBC), care must be taken with these
    for i=1:length(pos) % ONLY real atoms, no pad atoms, this is different in neighborlist 
        if abs(pos(i,1)-xbandl)<=10^-15
            band_atoms(end+1) = i;
        end
        
        if abs(pos(i,1)-xbandr)<=10^-15
            band_atoms(end+1) = i;
        end
    end
    band_atoms_mean = zeros(length(band_atoms),2);
    band_atoms_pos = zeros(size(band_atoms_mean));
    for i = 1:length(band_atoms)
        band_atoms_pos(i,:) = refpos(band_atoms(i),:);
        band_atoms_mean(i,:) = refpos(band_atoms(i),:);
    end
        
    %Associate a node to every band atom (needed as a reference for the new algorithm)
    disp('Associating Band Atoms with Interface Atoms ...')
    band_atom_ass_atoms = associateIntAtomswBandAtoms(example_name, band_atoms, pos, fe_model,  xbandl, xbandr, 0, 0, interface_atom_list);
    band_atom_ass_atoms_position = zeros(length(band_atom_ass_atoms),2);
  
    
    elseif strcmp(example_name, 'Tensile_Test_rotated')
    band_atoms = [];
    band_atoms_mean = [];
    lower_edge_atoms = [];
    pbc_atoms = [];
    
    xbandl = round(bandatom_position)*dx;
    xbandr = xbox - xbandl;
    
    pbc_atoms = []; %list of ''excluded'' atoms (cause of PBC), care must be taken with these
    for i=1:length(pos) % ONLY real atoms, no pad atoms, this is different in neighborlist 
        if abs(pos(i,1)-xbandl)<=10^-15
            band_atoms(end+1) = i;
        end
        
        if abs(pos(i,1)-xbandr)<=10^-15
            band_atoms(end+1) = i;
        end
    end
    band_atoms_mean = zeros(length(band_atoms),2);
    band_atoms_pos = zeros(size(band_atoms_mean));
    for i = 1:length(band_atoms)
        band_atoms_pos(i,:) = refpos(band_atoms(i),:);
        band_atoms_mean(i,:) = refpos(band_atoms(i),:);
    end
        
    %Associate a node to every band atom (needed as a reference for the new algorithm)
    disp('Associating Band Atoms with Interface Atoms ...')
    band_atom_ass_atoms = associateIntAtomswBandAtoms(example_name, band_atoms, pos, fe_model,  xbandl, xbandr, 0, 0, interface_atom_list);
    band_atom_ass_atoms_position = zeros(length(band_atom_ass_atoms),2);
    
end

i=0;
for j = band_atom_ass_atoms
    i = i+1;
    band_atom_ass_atoms_position(i,:) = pos(j,:);
end


if test_band_atom_association
    % TEST BAND POSITION AND ASSOCIATION
    figure(30)
    scatter(pos(band_atoms,1),pos(band_atoms,2))
    hold on
    scatter(band_atom_ass_atoms_position(:,1),band_atom_ass_atoms_position(:,2))
    axis equal
end


if strcmp(example_name, 'Example2')
    %set inactive atoms (Example2, top atoms and right atoms)
    refposall = [refpos; refpadpos];
    for i=1:length(refposall(:,1))
        if abs(refposall(i,2)-ybox)<=10^-15 %careful, comparison with ybox is only valid if ymin = 0!
            exclude_atoms(end+1) = i;
        end
        if refposall(i,1)>=axmax-5.1*dx
            inactive_atoms(end+1) = i;
        end
    end
    exclude_atoms = unique(exclude_atoms);
    inactive_atoms = unique(inactive_atoms);
    elements = [];
    virial_stress_atoms=[];
    for i=1:length(refpos(:,1))
        if isempty(find(inactive_atoms==i)) && isempty(find(exclude_atoms==i))
            if refpos(i,1) > axmin + 5*dx-tol
                virial_stress_atoms(end+1) = i;
            end
        end
    end
elseif strcmp(example_name, 'Dislocation')
    exclude_atoms = [];
    lower_edge_atoms = [];
    for i = 1:length(pos(:,1))
        %find top/bottom atoms
        if pos(i,2) >= max(pos(:,2))-1.01*rcut
            if pos(i,1) >= (axmax+axmin)/2-2*a && pos(i,1) <= (axmax+axmin)/2+2*a;
                inactive_atoms(end+1) = i;
                if pos(i,2) <= max(pos(:,2))-rcut + 0.01*a
                    lower_edge_atoms(end+1) = i;
                end
            end
        end
    end
    elements = CreateDetectionBand(12*dy,a);
    virial_stress_atoms = [];
    burgers_history = zeros(2, length(elements));
    burgers_time_record = zeros(1,length(elements));
elseif strcmp(example_name, 'Radial_Pulse')
    exclude_atoms = [];
    lower_edge_atoms = [];
    for i = 1:length(pos(:,1))
        %find top/bottom atoms
        midx = (axmax+axmin)/2;
        midy = (aymax+aymin)/2;
        
        
        if sqrt( (pos(i,1)-midx)^2 + (pos(i,2)-midy)^2 )<=rcut
            if pos(i,1) >= (axmax+axmin)/2-2*a && pos(i,1) <= (axmax+axmin)/2+2*a
                inactive_atoms(end+1) = i;
            end
        end
        
    end
    elements = [];
    virial_stress_atoms = [];
    
elseif strcmp(example_name, 'Example2_rotated')
    %set inactive atoms (Example2, top atoms and right atoms)
    refposall = [refpos; refpadpos];
    virial_stress_atoms = [];
    for i=1:length(refposall(:,1))
        if abs(refposall(i,2)-ybox)<=10^-15 %careful, comparison with ybox is only valid if ymin = 0!
            exclude_atoms(end+1) = i;
        end
        if refposall(i,1)>=axmax-5.1*dx
            inactive_atoms(end+1) = i;
        end
    end
    exclude_atoms = unique(exclude_atoms);
    inactive_atoms = unique(inactive_atoms);
    elements = [];
    for i=1:length(refpos(:,1))
        if isempty(find(inactive_atoms==i)) && isempty(find(exclude_atoms==i))
            if refpos(i,1) > axmin + 5*dx-tol
                virial_stress_atoms(end+1) = i;
            end
        end
    end
    
elseif strcmp(example_name, 'Tensile_Test')
    %set inactive atoms (Example2, top atoms and right atoms)
    exclude_atoms = [];
    inactive_atoms = [];
    elements = [];
    virial_stress_atoms = [];
    
    %alle Atome
    for i=1:length(refpos(:,1))
	if isempty(find(inactive_atoms==i)) && isempty(find(exclude_atoms==i))
        if refpos(i,1) > axmin + 5*a && refpos(i,1) < axmax - 5*a
            virial_stress_atoms(end+1) = i;
        end
	end
    end
    
    elseif strcmp(example_name, 'Tensile_Test_rotated')
    %set inactive atoms (Example2, top atoms and right atoms)
    exclude_atoms = [];
    inactive_atoms = [];
    elements = [];
    virial_stress_atoms = [];
    
    %alle Atome
    for i=1:length(refpos(:,1))
	if isempty(find(inactive_atoms==i)) && isempty(find(exclude_atoms==i))
        if refpos(i,1) > axmin + 15*dx-tol && refpos(i,1) < axmax - 15*dx+tol
            if refpos(i,2) > 6*dy+rcut*1.1 && refpos(i,2) < 11*dy-rcut*1.1
                virial_stress_atoms(end+1) = i;
            end
        end
	end
    end

end


end

