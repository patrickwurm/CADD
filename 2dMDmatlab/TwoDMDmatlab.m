function TwoDMDmatlab(fe_model, fe_analysis, settings, sim_name)

close all;
example_name = fe_analysis.example;

if isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Analysis')
    cont_type = 'static';
elseif isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Dynamic_Analysis')
    cont_type = 'dynamic';
elseif isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Hybrid_Analysis')
    cont_type = 'hybrid';
end

rng('default')
if strcmp(settings.additional_input.example_id_1,'Tensile_Test')
    rng(settings.additional_input.rng_seed);
else
    rng(1);
end

% UNIT SYSTEM (SI)
% (EAM Parameters must be converted)
% mass: kg
% distance: m
% time: s
% energy: J
% velocity: m/s
% temperature: Kelvin

if isfield(settings.additional_input,'temperature')
    if settings.additional_input.temperature == 0.001
        a = 2.6485095*1e-10;
    elseif settings.additional_input.temperature == 50
        a = 2.6491000*1e-10;
    elseif settings.additional_input.temperature == 100
        a = 2.6497690*1e-10; %OLD
        %a = 2.6492*1e-10; %NEW
    elseif settings.additional_input.temperature == 150
        a = 2.6507000*1e-10;
    elseif settings.additional_input.temperature == 250
        a = 2.6530000*1e-10;
    else
        error('No material parameters saved for target temperature.')
    end
else
    a = 2.6485095*1e-10; %0.001K EAM
end

fe_model.setA(a);
fe_model.bc_handler.incorporateBC(fe_model, 0, 0); %this is just a first call to initiate some quantities in the bc handler

rcut = 6.403928746393108*1e-10;
axmin = 0;
aymin = 0;

if strcmp(example_name,'Example2_rotated')
    dx = a;
    dy = sqrt(3)*a;
    if strcmp(settings.additional_input.example_id_1,'Equilibrium_small')
        axmax = 100*dx;
    else
        axmax = 500*dx;
    end
    aymax = 3*dy;
elseif strcmp(example_name,'Example2')
    dx = sqrt(3)*a;
    dy = a;
    axmax = 289*dx;
    aymax = 6*dy;
elseif strcmp(example_name,'Tensile_Test')
    dx = sqrt(3)*a;
    dy = a;
    axmax = 62*dx;
    aymax = 29*dy;
elseif strcmp(example_name,'Tensile_Test_rotated')
    dx = a;
    dy = sqrt(3)*a;
    axmax = 110*dx;
    aymax = 17*dy;
end

global tol
tol = 1e-16;

global refpos;
global refpadpos;
global inactive_atoms;
global exclude_atoms;
global virial_stress_atoms;
exclude_atoms = [];
inactive_atoms = [];

kb = 1.3806485279e-23; %J/K
amass = 26.9815385*1.66053904020e-27; %kg

settings.continuum_external_change = [1,1];
settings.continuum_external_displacement = 0*a;
settings.write_reflection_file = 0;
restart_sim = 0;
mstep_save_at = 0;
if strcmp(settings.additional_input.example_id_1,'Reflection')
    if settings.additional_input.pulse == 0
        mstep = 20100;
        mstep_save_at = 20000;
        restart_sim = 0;
        settings.write_reflection_file = 0;
    elseif settings.additional_input.pulse == 1
        mstep = 39000;
        mstep_save_at = 0;
        restart_sim = 1;
        settings.int_forth = [1, 20000];
        settings.int_back = [34000, mstep-100];
        settings.write_reflection_file = 1;
    elseif settings.additional_input.pulse == 2
        mstep = 50000;
        mstep_save_at = 0;
        restart_sim = 1;
        settings.int_forth = [1, 29000];
        settings.int_back = [45000, mstep-100];
        settings.write_reflection_file = 1;
    elseif settings.additional_input.pulse == 3
        mstep = 55000;
        mstep_save_at = 0;
        restart_sim = 1;
        settings.int_forth = [1, 29000];
        settings.int_back = [50000, mstep-100];
        settings.write_reflection_file = 1;
    else
        error('Unknown Pulse')
    end
elseif strcmp(settings.additional_input.example_id_1,'Tensile_Test')
    if settings.additional_input.temperature == 50
            mstep = 160000;
            settings.continuum_external_change = [40001,160000];
            settings.continuum_external_displacement = 8.58*a;
    elseif settings.additional_input.temperature == 150
            mstep = 160000;
            settings.continuum_external_change = [40001,160000];
            settings.continuum_external_displacement = 8.58*a;
    elseif settings.additional_input.temperature == 250
            mstep = 200000;
            settings.continuum_external_change = [40001,200000];
            settings.continuum_external_displacement = 11.44*a;
    end
elseif strcmp(settings.additional_input.example_id_1,'Equilibrium')
    if settings.additional_input.speed == 1
        mstep = 400000;
        settings.continuum_external_change = [20001,20001+36363];
        settings.continuum_external_displacement = -10*a;
    elseif settings.additional_input.speed == 2
        mstep = 400000;
        settings.continuum_external_change = [20001,20001+12121];
        settings.continuum_external_displacement = -10*a;
    elseif settings.additional_input.speed == 3
        mstep = 400000;
        settings.continuum_external_change = [20001,20001+4040];
        settings.continuum_external_displacement = -10*a;
    end
elseif strcmp(settings.additional_input.example_id_1,'Equilibrium_small')
    if settings.additional_input.speed == 1
        mstep = 50000;
        settings.continuum_external_change = [20001,20001+12121];
        settings.continuum_external_displacement = -2*a;
    end
end

%mstep = 60000; %total number of time steps
mstep_run_in = 10^6; %time period to record the equilibrium fluctutions of the demand-based algorithm
%mstep_save_at = 0; %optional: time step at which the current state should be saved for later runs (restart_sim)
%restart_sim = 1; %1: load initial state from file (saved before at mstep_save_at); 0: run simulation without loading anything
load_pad_association = 0; %load pad atom to element association from file (to save time)
load_fluc = 0; %load existing fluctuations from file

iplot = 0; %plot atomic positions, energies, temperature, .. every iplot steps
inl = 100; %update neighborlist every inl steps
iwriteVTK = 50; %write configuration to VTK file every iwriteVTK steps
irunDD = 100; %update positions of continuum dislocations every irunDD timesteps
fe_step = 0; %this counts how often the FEA was invoked (this is for output purposes)

if isfield(settings.additional_input,'temperature')
    tset = settings.additional_input.temperature; %target temperature in K
else
    tset = 0.001; %target temperature in K
end

%Initialize NH thermostat
xi=1; %xi, xi1, xi2: initial values for xi in the NH thermostat
xi1=0;
xi2=0;
tau = 5e-14;

deltat = fe_model.time_step;

if strcmp(example_name,'Example2') || strcmp(example_name,'Example2_rotated')
    pbc = true; %this is tricky and depends on the example, be careful
else
    pbc = false;
end
D = 40; %downsampling factor used in the demand-based algorithm

M_window = round(1500/D); %Window size
fc = 1.885e12/(2*pi)/(1e15/D); %cutoff frequency

if strcmp(settings.additional_input.example_id_1,'Reflection')
    settings.thermostatted = 0; %0: thermostat off, 1: NH, 2: stadium damping/SBC 3: SBC (no stadium)
else
    settings.thermostatted = 3; %0: thermostat off, 1: NH, 2: stadium damping/SBC 3: SBC (no stadium)
end

settings.atomic_potential = 4; %1: EAM, 2: LJ, 3:LJ (vectorized), 4:EAM (vectorized)
settings.test_band_atom_association = 0;
settings.print_detection_band = 0;
settings.include_continuum_dislocations = 0;
settings.recover_velocity_center_of_mass = 0; %be very careful with this one
settings.plot_positions = 0;
settings.save_workspace_after_finishing = 1;
settings.damping_band_depth = 2.1*dx;
settings.gamma_0 = 2.5*10^13;
settings.effective_cutoff = 1.10*rcut; %effective cutoff radius used for the neighborlists

settings.output_handler = fe_analysis.output_handler;

dislocation_grow_up_time = 500;

%initialize Rainbow Tables
if settings.atomic_potential == 4
    force_tables = initializeRainbowTables(a, rcut,settings.atomic_potential);
else
    force_tables = [];
end

disl_drag_coeff = 0.75*1e-4;

%reference vectors used in the dislocation detection
reference_vectors = [ 0 ,sqrt(3)/2*a, sqrt(3)/2*a, 0, -sqrt(3)/2*a, -sqrt(3)/2*a, -sqrt(3)/6*a ,sqrt(3)/6*a, sqrt(3)/3*a, sqrt(3)/6*a, -sqrt(3)/6*a, -sqrt(3)/3*a; ...
    a, 1/2*a, -1/2*a, -a, -1/2*a, 1/2*a, 1/2*a, 1/2*a, 0, -1/2*a -1/2*a, 0]; %

delta = 20*a; %transition distance of continuum_dislocations

if restart_sim == 1
    mstep_save_at = 0;
end
if load_fluc == 1
    mstep_run_in = 0;
end

restart_settings.istep=1;
restart_settings.mstep=mstep;
restart_settings.iwriteVTK = iwriteVTK;
restart_settings.mstep_save_at = mstep_save_at;
restart_settings.mstep_run_in = mstep_run_in;
restart_settings.restart_sim = 0;
restart_settings.i_percent = 1;
restart_settings.inv_cont_for = 0;
restart_settings.load_fluc = load_fluc;
restart_settings.trigger_case = settings.trigger_case;
restart_settings.settings = settings;

% =========================================================================

disp(' ')
disp('I am creating the MD Model ...')
if strcmp(example_name,'Example2_rotated') || strcmp(example_name,'Tensile_Test_rotated')
    pos = genatoms_original(a, axmin, aymin, axmax, aymax, dx, dy); %Generate lattice
else 
    pos = genatoms(a, axmin, aymin, axmax, aymax, 0, dx, dy); %Generate lattice
end

padpos = genpadatoms(a, axmin, aymin, axmax, aymax, rcut, example_name, dx, dy);

if 2>1
    if strcmp(example_name,'Tensile_Test') || strcmp(example_name,'Tensile_Test_rotated')
        [pos, padpos] = delete_atoms(a, axmin, aymin, axmax, aymax, pos, padpos, example_name, dx, dy);
    end
end

mxatm = length(pos(:,1));
refpos = pos;

mxpatm = length(padpos(:,1));
refpadpos = padpos;

interface_atom_list = findInterfaceAtoms(refpos, axmin, aymin, axmax, aymax, a, example_name);
interface_node_list = associateIntAtomswNodes(fe_model, refpos, interface_atom_list);

if strcmp(example_name,'Example1') || strcmp(example_name,'Dislocation') || strcmp(example_name,'Radial_Pulse') % is this right for the dislocation example?
    boundary_cond_matrix_interface = zeros(fe_model.dimension*length(interface_node_list),6);
elseif strcmp(example_name,'Example2') || strcmp(example_name,'Example2_rotated') || strcmp(example_name,'Tensile_Test') || strcmp(example_name,'Tensile_Test_rotated')
    boundary_cond_matrix_interface = zeros(length(interface_node_list),6);
else
    error('Example not found')
end

if load_pad_association
    disp('Loading pad association ...')
    disp('-> Caution: If this is the first run of the example, you need to create and save a pad association first to load for later runs. Set load_pad_association to 0 and restart.')
    load(['examples/',fe_analysis.example,'/pad_assoc/pad_assoc.mat'],'fe_model')
    fe_analysis.setModel(fe_model);
    disp('... Done loading')
else
    disp('Associating Pad Atoms with Elements ...')
    associatePadAtomswElements(fe_model, refpadpos, mxatm);
    %disp('Saving pad associations for later runs ...')
    %mkdir(['examples/',fe_analysis.example,'/pad_assoc'])
    %save(['examples/',fe_analysis.example,'/pad_assoc/pad_assoc.mat'],'fe_model','-v7.3')
    %disp('... Done saving')
end

findNatCoordsPad(fe_model);

fe_analysis.ComputeHMatrix(mxpatm, mxatm)

xbox = max(pos(:,1))-min(pos(:,1)); %size of the atomistic domain in x-direction
ybox = max(pos(:,2))-min(pos(:,2)); %size of the atomistic domain in y-direction

if strcmp(example_name,'Tensile_Test') || strcmp(example_name,'Tensile_Test_rotated') 
    interface_atom_list_left = [];
    interface_atom_list_right = [];
    for i=1:length(interface_atom_list)
        if (abs(pos(interface_atom_list(i),1)-0)<tol) 
		interface_atom_list_left(end+1) = interface_atom_list(i);
	elseif (abs(pos(interface_atom_list(i),1)-xbox)<tol)
		interface_atom_list_right(end+1) = interface_atom_list(i);
        end
    end
end


[pbc_atoms, band_atoms, band_atoms_pos, band_atom_ass_atoms, band_atom_ass_atoms_position, inactive_atoms, exclude_atoms, lower_edge_atoms, elements, virial_stress_atoms] = getSpecialAtoms(example_name, pos, refpos, refpadpos, xbox, ybox, settings.bandatom_position, a, axmin, axmax, aymin, aymax, rcut, inactive_atoms, exclude_atoms, fe_model, interface_atom_list, settings.test_band_atom_association, dx, dy);

if settings.trigger_case == 7
    [kernel_matrix, pos_matrix] = initializeConvolutionKernel(M_window, band_atoms, band_atom_ass_atoms, fc, D, pos );
end

run_in.rel_distance = zeros(length(band_atom_ass_atoms),2,floor(mstep_run_in/D));
run_in.min_distance = zeros(length(band_atom_ass_atoms),2);
run_in.max_distance = zeros(length(band_atom_ass_atoms),2);

band_atoms_pos_hist = zeros(length(band_atoms),2,mstep+1);
band_atoms_pos_hist(:,:,1) = band_atoms_pos(:,:);
total_dislacement_hist = zeros(1, mstep);
max_disp_x_hist = zeros(1, mstep);

edge_force_hist = zeros(1,mstep+1);
edge_force_hist(1) = 0;
    
interface_position_hist = zeros(1,mstep);
strain_hist = zeros(1,mstep);
strain_hist2 = zeros(1,mstep);

mean_temp_container = zeros(1,mstep);
temp_container = zeros(1,mstep);
temp2_container = zeros(1,mstep);
temp_inner_container = zeros(1,mstep);
temp_outer_container = zeros(1,mstep);

cont_on = zeros(1,mstep);
cont_off = zeros(1,mstep);
inv_cont_for = 0;

vel = zeros(mxatm,2);
force = zeros(mxatm,2);

vel = intializeMD(vel, tset, kb, amass, mxatm, refpos, aymax, a);

if iplot > 0
figure(1)
plotatoms = scatter(pos(:,1),pos(:,2),'o','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'SizeData',30);
hold on 
plotpadatoms = scatter(padpos(:,1),padpos(:,2),'o','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 0 0],'SizeData',30);
allpos = [pos; padpos];
plotdisloc(1) = scatter([],[]); %placeholder
hold on
xlim([min(pos(:,1))-20*a max(pos(:,1))+20*a])
ylim([min(pos(:,2))-60*a max(pos(:,2))+2*a])
axis equal


figure(2)
subplot(3,1,1)
plotenergy = plot(0,0);
grid on
title('Energy vs. time')
subplot(3,1,2)
plottemp = plot(0,0);
grid on
title('Instant. temperature vs. time')
subplot(3,1,3)
plotmeantemp = plot(0,0);
grid on
title('Mean temperature vs. time')
end

disp('Writing "standard" vtk output file ...')
vtks = 0;
fe_analysis.output_handler.write_standard(fe_model, vtks, [refpos; refpadpos], [pos; padpos], mxatm, mxpatm, 0);

%output referential configuration
if ~isempty(fe_analysis.output_handler)
    energies = zeros(length(pos(:,1))+mxpatm,1);
    fe_analysis.output_handler.write(fe_model, vtks, [refpos; refpadpos], [pos; padpos], mxatm, mxpatm, 0, energies, vel, fe_analysis.u1);
    vtks = vtks + 1;
end

for i = 1:length(exclude_atoms)
    force(exclude_atoms(i),:) = zeros(1,2);
    vel(exclude_atoms(i),:) = zeros(1,2);
end

for i = 1:length(inactive_atoms)
    force(inactive_atoms(i),:) = zeros(1,2);
    vel(inactive_atoms(i),:) = zeros(1,2);
end

dump = initializePosDump(example_name, mstep, refpos, fe_model);


if strcmp(example_name,'Tensile_Test') || strcmp(example_name,'Tensile_Test_rotated') || strcmp(example_name,'Example2_rotated')
    dump.mean_stress_xx = zeros(1,mstep);
    dump.mean_stress_yy = zeros(1,mstep);
    dump.mean_stress_xy = zeros(1,mstep);
    dump.RMSD = zeros(1,mstep);
    dump.RMSD_ref = zeros(1,mstep);
    dump.AAD = zeros(1,mstep);
    dump.AAD_ref = zeros(1,mstep);
end

%Update Neighborlist
[ nl, numneigh ] = neighborlist([pos; padpos], mxatm, mxpatm, rcut, pbc, xbox, ybox, a, example_name, settings);

if restart_sim == 1
    close all
    example = fe_analysis.example;
    clearvars -except restart_settings example cont_type tset
    restart_settings0 = restart_settings;
    disp('===================================')
    disp('Continuing simulation from file (restart_sim == 1) ...')
    load(['examples/',example,'/restart_sim/workspace_',cont_type,'_',num2str(tset),'K.mat'])
    restart_settings = restart_settings0;
    clear restart_settings0 example
    [istep, mstep, iwriteVTK, mstep_save_at, mstep_run_in, restart_sim, i_percent, timer_1, inv_cont_for, load_fluc, settings] = initializeRestart(fe_model, fe_analysis, restart_settings);
    interface_position_hist = zeros(1,mstep);
    strain_hist = zeros(1,mstep);
    strain_hist2 = zeros(1,mstep);
    fe_analysis.output_handler = settings.output_handler;
    clear band_atoms_pos_hist
    fe_analysis.output_handler = settings.output_handler;
    
    if strcmp(settings.additional_input.example_id_1,'Reflection')
        fe_analysis.u1 = zeros(size(fe_analysis.u1));
    end
end

atom_load = setupLoad(example_name, a, axmin, axmax, aymin, aymax, rcut, settings);

if ~isfield(settings.additional_input,'fe_damping_coefficient')
    settings.additional_input.fe_damping_coefficient=0;
end

if isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Dynamic_Analysis')
    if isfield(settings.additional_input,'selective_fe_damping')
        if settings.additional_input.selective_fe_damping
            damp_node_indices_del = [];
            damp_matrix_indices_del = [];
            for node = fe_model.node_dict
                if node.material_coordinates(1) < -6*a
                    damp_node_indices_del(end+1) = node.number;
                end
            end
            for i = 1:length(damp_node_indices_del)
                damp_matrix_indices_del(end+1) = (damp_node_indices_del(i)*2)-1;
                damp_matrix_indices_del(end+1) = (damp_node_indices_del(i)*2);
            end
            fe_analysis.damp_matrix = fe_analysis.mass_matrix * settings.additional_input.fe_damping_coefficient;
            fe_analysis.damp_matrix(damp_matrix_indices_del,:) = 0;
            fe_analysis.damp_matrix = sparse(fe_analysis.damp_matrix);
        else
            fe_analysis.damp_matrix = sparse(fe_analysis.mass_matrix * settings.additional_input.fe_damping_coefficient);
        end
    end
    
elseif isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Hybrid_Analysis')
    if isfield(settings.additional_input,'selective_fe_damping')
        if settings.additional_input.selective_fe_damping
            damp_node_indices_del = [];
            damp_matrix_indices_del = [];
            for node = fe_model.node_dict
                if node.material_coordinates(1) < -6*a
                    damp_node_indices_del(end+1) = node.number;
                end
            end
            for i = 1:length(damp_node_indices_del)
                damp_matrix_indices_del(end+1) = (damp_node_indices_del(i)*2)-1;
                damp_matrix_indices_del(end+1) = (damp_node_indices_del(i)*2);
            end
            fe_analysis.dynamic_analysis.damp_matrix = fe_analysis.dynamic_analysis.mass_matrix * settings.additional_input.fe_damping_coefficient;
            fe_analysis.dynamic_analysis.damp_matrix(damp_matrix_indices_del,:) = 0;
            fe_analysis.dynamic_analysis.damp_matrix = sparse(fe_analysis.dynamic_analysis.damp_matrix);
        else
            fe_analysis.dynamic_analysis.damp_matrix = sparse(fe_analysis.dynamic_analysis.mass_matrix * settings.additional_input.fe_damping_coefficient);
        end
    end
end
    
if strcmp(example_name,'Tensile_Test') || strcmp(example_name,'Tensile_Test_rotated')
    npart = 20;
    area_split = initializeAreaSplit(refpos, npart, a, xbox, ybox, example_name, dx, dy);
end

if isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Analysis')
    xtsample = zeros(settings.delta_t_static,length(interface_node_list));
    ytsample = zeros(settings.delta_t_static,length(interface_node_list));
    avgcount = 0;
end

disp(' ')
disp('Running Simulation ...')
disp(' ')
avg_counter=0;
i_percent = 1;
timer_0 = tic;
timer_1 = tic;
updateDD = 0;
write_to_nodes = 0;
istep = 0;
cont_is_off = true;
last_off = 0;

if strcmp(settings.additional_input.example_id_1,'Equilibrium') || strcmp(settings.additional_input.example_id_1,'Equilibrium_small')
    if strcmp(settings.additional_input.example_id_1,'Equilibrium')
        ref_indices = [];
        for i = 1:mxatm
            if refpos(i,1) < 494.99*a
                ref_indices(end+1) = i;
            else
            end
        end
        u_desired = 10/995*(495*a - pos(ref_indices,1));
    elseif strcmp(settings.additional_input.example_id_1,'Equilibrium_small')
        ref_indices = [];
        for i = 1:mxatm
            if refpos(i,1) < 94.99*a
                ref_indices(end+1) = i;
            else
            end
        end
        u_desired = 2/195*(95*a - pos(ref_indices,1));
    end
end

while istep < mstep

    if strcmp(settings.additional_input.example_id_1,'Equilibrium') || strcmp(settings.additional_input.example_id_1,'Equilibrium_small')
            if istep < 2500
                settings.gamma_0 = 2.5*10^13;
                settings.damping_band_depth = 475*dx;
            elseif istep < 5000
                settings.gamma_0 = 2.5*10^13;
                settings.damping_band_depth = 25*dx;
            else
                settings.damping_band_depth = 2.1*dx;
                settings.gamma_0 = 3/4*2.5*10^13;
            end
    end
    
    if strcmp(settings.additional_input.example_id_1,'Tensile_Test') 
        if istep < 2500
            settings.gamma_0 = 2.5*10^13;
            settings.damping_band_depth = 50*dx;
        elseif istep < 5000
            settings.gamma_0 = 2.5*10^13;
            settings.damping_band_depth = 25*dx;
        else
            settings.damping_band_depth = 2.1*dx;
            settings.gamma_0 = 2.5*10^13;
            %settings.gamma_0 = 4.52*10^12;
        end
    end

    istep = istep+1;    
    avg_counter = avg_counter+1;
    
 
    if istep==1
        damped_atoms = [];
        for ia = 1:mxatm
            if refpos(ia,1) <= axmin + settings.damping_band_depth
                damped_atoms(end+1) = ia;
            end
        end
    end

    fe_model.setTime(istep);
    
    if abs(istep/(round(mstep/100)*100)*100-i_percent*1)<=10^-6
        
        %disp([num2str(istep/mstep*100),'% (', num2str(istep),' steps) finished in ',num2str(toc(timer_1)),' seconds.'])
        
        estimated_time_left = toc(timer_0)*(mstep/istep-1)/60;
        term_out = sprintf('%5.1f%% (%i steps) finished in %5.2f seconds. Estimated time left: %6.1f minutes.',istep/mstep*100,istep,toc(timer_1),estimated_time_left);
        disp(term_out)
        timer_1 = tic;
        i_percent=i_percent+1;
    end
    
    pos = applyLoad(inactive_atoms, atom_load, pos, refpos, example_name, istep, settings);
    
    [pos, vel, xi, xi1, xi2, force, ekin, epot, temp, energies, virial] = Solve(pos, padpos, vel, xi, xi1, xi2, force, deltat, tau, mxatm, mxpatm, amass, rcut, xbox, ybox, tset, pbc, kb, nl, numneigh, a, settings, example_name, force_tables);
    
    interface_position_hist(istep) = mean(pos(interface_atom_list,1));
    
    if strcmp(example_name, 'Tensile_Test') || strcmp(example_name, 'Tensile_Test_rotated')
        strain_hist(istep) = (mean(pos(interface_atom_list_right,1))-mean(pos(interface_atom_list_left,1))-xbox)/xbox;
    end
    
    [total_disp, max_disp_x] = calcTotalDisplacement(pos,refpos);
    
    total_dislacement_hist(istep) = total_disp/mxatm/a;
    max_disp_x_hist(istep) = max_disp_x/a;
        
    if strcmp(settings.additional_input.example_id_1,'Equilibrium') || strcmp(settings.additional_input.example_id_1,'Equilibrium_small')
        u_vec = pos(ref_indices,1)- refpos(ref_indices,1);
        u_dev = abs(u_vec - u_desired);
        u_dev2 = (u_vec - u_desired).^2;
        dump.AAD(istep) = 1/length(u_dev)*sum(u_dev);
        dump.RMSD(istep) = sqrt(1/length(u_dev)*sum(u_dev2));
        if istep < 20000
            dump.AAD_ref(istep) = 1/length(u_vec)*sum(abs(u_vec));
            dump.RMSD_ref(istep) = sqrt(1/length(u_vec)*sum(u_vec.^2));
        end
    end
    
    if strcmp(example_name, 'Example2') || strcmp(example_name, 'Example2_rotated')
        %project positions, velocities and forces of the real atoms to the corresponding pbc_atoms
        for i=1:length(pbc_atoms(1,:)) % ONLY real atoms, no pad atoms, this is different in neighborlist
            pos(pbc_atoms(1,i),:) = refpos(pbc_atoms(1,i),:) + (pos(pbc_atoms(2,i),:) - refpos(pbc_atoms(2,i),:));
            vel(pbc_atoms(1,i),:) = vel(pbc_atoms(2,i));
            force(pbc_atoms(1,i),:) = force(pbc_atoms(2,i),:);
        end
        
        v1m = mean(vel(virial_stress_atoms,1));
        v2m = mean(vel(virial_stress_atoms,2));
        %v1m = 0; %Kein Abzug von Mittelwert in per atom Definition? Stimmt das?
        %v2m = 0; %Kein Abzug von Mittelwert in per atom Definition? Stimmt das?
        if settings.atomic_potential == 4
            for i=1:length(virial_stress_atoms)
                ia = virial_stress_atoms(i);
                virial.spa_xx(i) = -amass*(vel(ia,1)-v1m)*(vel(ia,1)-v1m) - virial.sumvirial_xx(ia);
                virial.spa_yy(i) = -amass*(vel(ia,2)-v2m)*(vel(ia,2)-v2m) - virial.sumvirial_yy(ia);
                virial.spa_xy(i) = -amass*(vel(ia,1)-v1m)*(vel(ia,2)-v2m) - virial.sumvirial_xy(ia);
            end
        else
            for i=1:length(virial_stress_atoms)
                ia = virial_stress_atoms(i);
                virial.spa_xx(i) = -amass*(vel(ia,1)-v1m)*(vel(ia,1)-v1m) - virial.sumvirial_xx(i);
                virial.spa_yy(i) = -amass*(vel(ia,2)-v2m)*(vel(ia,2)-v2m) - virial.sumvirial_yy(i);
                virial.spa_xy(i) = -amass*(vel(ia,1)-v1m)*(vel(ia,2)-v2m) - virial.sumvirial_xy(i);
            end
        end
        
        if istep==1 || ~mod(istep,10)
            [k, stress_area] = boundary(pos(virial_stress_atoms,1),pos(virial_stress_atoms,2));
        end
        
        dump.mean_stress_xx(istep) = 1/stress_area*sum(virial.spa_xx);
        dump.mean_stress_yy(istep) = 1/stress_area*sum(virial.spa_yy);
        dump.mean_stress_xy(istep) = 1/stress_area*sum(virial.spa_xy);
        
    elseif strcmp(example_name, 'Dislocation')
        pbc_atoms = [];
        edge_force_hist(istep+1) = mean(force(lower_edge_atoms,2));
    elseif strcmp(example_name, 'Tensile_Test') || strcmp(example_name, 'Tensile_Test_rotated') 
        v1m = mean(vel(virial_stress_atoms,1));
        v2m = mean(vel(virial_stress_atoms,2));
        for i=1:length(virial_stress_atoms)
            ia = virial_stress_atoms(i);
            virial.spa_xx(i) = -amass*(vel(ia,1)-v1m)*(vel(ia,1)-v1m) - virial.sumvirial_xx(ia);
            virial.spa_yy(i) = -amass*(vel(ia,2)-v2m)*(vel(ia,2)-v2m) - virial.sumvirial_yy(ia);
            virial.spa_xy(i) = -amass*(vel(ia,1)-v1m)*(vel(ia,2)-v2m) - virial.sumvirial_xy(ia);
        end
        if istep==1 || ~mod(istep,10)
            [k, stress_area] = boundary(pos(virial_stress_atoms,1),pos(virial_stress_atoms,2));
        end
        
        [lowest, indxlowest] = mink(pos(virial_stress_atoms,1),6);
        [highest, indxhighest] = maxk(pos(virial_stress_atoms,1),6);
        
        ml = mean(lowest);
        mh = mean(highest);
        
        if istep==1
            ml0 = ml;
            mh0 = mh;
            l0 = mh0 - ml0;
        end
        
        l1 = mh - ml;
        strain_hist2(istep) = (l1-l0)/l0;
        
        dump.mean_stress_xx(istep) = 1/stress_area*sum(virial.spa_xx);
        dump.mean_stress_yy(istep) = 1/stress_area*sum(virial.spa_yy);
        dump.mean_stress_xy(istep) = 1/stress_area*sum(virial.spa_xy);
    else
        %do nothing
    end
    
    if settings.trigger_case == 7
        [inv_cont_for, trigger, run_in, pos_matrix, band_atom_ass_atoms_position, filter_output] = AlgorithmControl(settings.trigger_case, pos_matrix, kernel_matrix, band_atoms, band_atom_ass_atoms, band_atom_ass_atoms_position, run_in, mstep_run_in, fe_analysis.example, load_fluc, restart_sim, D, istep, inv_cont_for, pos);
    else
        [trigger, filter_output] = Cont_trigger_DSP(settings.trigger_case);
    end

    additional_displacements_pad = zeros(length(padpos(:,1)),2); %these are the u_tilde displacements on the pad atoms

    %Be careful, this part is highly specific and has not been updated for
    %the hybrid continuum
    if settings.include_continuum_dislocations
        if ~exist('dislocation_history','var')
            dislocation_history = zeros(length(elements(:,1)),1);
        end
        
        delta = 15*a; %Shifting distance of dislocation
        
        %CHECK IF ANY NEW DISLOCATIONS HAVE BEEN FOUND
        [ detected_in_el, burgers ] = DetectDislocations( reference_vectors, pos, elements, a, correction);
        
        %CREATE A NEW DIPOL IF THE FOUND DISLOCATIONS ARE NEW
        if 2>1
            if ~isempty(detected_in_el)
                detected_in_el
                burgers/a
                
                if ~mod(istep,10)
                    for i = 1:length(detected_in_el)
                        disp([num2str(istep),' ',num2str(detected_in_el(i)) ,' ', num2str(burgers(1,i)/a),' ',num2str(burgers(2,i)/a)])
                        elem_nums = reshape(elements(detected_in_el(i),:),[],1);
                        figure(77)
                        axis equal
                        for o = 1:length(elem_nums)
                            scatter(refpos(elem_nums(o),1),refpos(elem_nums(o),2), 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'SizeData',30 )
                            hold on
                            scatter(pos(elem_nums(o),1),pos(elem_nums(o),2), 'SizeData',30 )
                        end
                        drawnow
                    end
                end
                for i = 1:length(detected_in_el)

                    last_detected = burgers_time_record(detected_in_el(i));
                    
                    if istep >= last_detected + 2*dislocation_grow_up_time
                    else
                        continue
                    end

                    burgers(:,i) = [0; a*sign(burgers(2,i))];
                    burgers_time_record(detected_in_el(i)) = istep;
                    
                    disp(['Step: ',num2str(istep)])
                    disp(['Dislocation with burgers vector [',num2str(burgers(1,i)),',',num2str(burgers(2,i)),'] in element ',num2str(detected_in_el(i)),' detected.'])
                    
                    if detected_in_el(i) == 5
                        correction(detected_in_el(i),1) = 1;
                        correction(detected_in_el(i),2:3) = correction(detected_in_el(i),2:3) + [abs(burgers(1,i)), abs(burgers(2,i))];
                    elseif detected_in_el(i) == 8
                        correction(detected_in_el(i),1) = 2;
                        correction(detected_in_el(i),2:3) = correction(detected_in_el(i),2:3) + [abs(burgers(1,i)), abs(burgers(2,i))];
                    end
                    
                    
                    % Assumption: Dislocation is placed in the middle of the
                    % element
                    dislocation_coords = mean(refpos(elements(detected_in_el(i),:),:));

                    % Create new continuum dislocation
                    new_disl_number = length(fe_model.dislocation_dict)+1;
                    fe_model.addContinuumDislocation(nsModel.ContinuumDislocation(new_disl_number, dislocation_coords - [0, delta], sign(burgers(2,i)), fe_model.element_dict, istep, dislocation_grow_up_time, disl_drag_coeff, irunDD*deltat, irunDD));
                    % Create new ghost dislocation
                    ghost_number = length(fe_model.ghostdislocation_dict)+1;
                    fe_model.addGhostDislocation(nsModel.GhostDislocation(ghost_number, dislocation_coords, -sign(burgers(2,i)), istep, dislocation_grow_up_time));
                    
                    dislocation_history(detected_in_el(i)) = dislocation_history(detected_in_el(i))+1;
                    
                    write_to_nodes = 1;
                    trigger = 1;
                    iwriteVTK = 100;
                end
            end
        end
        
        disp_field_changed = 0;
        
        for i=1:length(fe_model.dislocation_dict)
            if istep < fe_model.dislocation_dict(i).birth_time + fe_model.dislocation_dict(i).grow_up_time
                disp_field_changed = 1;
                %Superpose dipol displacements (GHOST + NEW) to atomic displacements
                for j = 1:length(pos(:,1)) %this could be vectorized, but i dont think it is worth the time
                    pos(j,:) = pos(j,:) + fe_model.ghostdislocation_dict(i).getGhostIncrementalDisplacementAtAtom(pos(j,:), fe_model.a, fe_model.material_dict{1,1});
                    pos(j,:) = pos(j,:) + fe_model.dislocation_dict(i).getIncrementalDisplacementAtAtom(pos(j,:), fe_model.a, fe_model.material_dict{1,1});
                end
            end
            if istep == fe_model.dislocation_dict(i).birth_time + fe_model.dislocation_dict(i).grow_up_time
                %just finished, set finish flag (and associate pad atoms with
                %elements again
                disp(['Step: ',num2str(istep)])
                disp(['Finished passing of dislocation ',num2str(i)])
            end
            if 2>1
                %Superpose dipol displacements (GHOST + NEW) to pad atom displacements
                %CARE, THIS MUST ALWAYS BE DONE (ALSO AFTER PASSING, with full INTENSITY)
                fe_model.ghostdislocation_dict(i).computeIntensity(istep);
                fe_model.dislocation_dict(i).computeIntensity(istep);
                for j = 1:length(padpos(:,1)) %this could be vectorized, but i dont think it is worth the time
                    pad_disp = fe_model.ghostdislocation_dict(i).getGhostDisplacementAtAtom(refpadpos(j,:), fe_model.a, fe_model.material_dict{1,1});
                    %padpos(j,:) = padpos(j,:) + pad_disp;
                    additional_displacements_pad(j,:) = additional_displacements_pad(j,:) + pad_disp;
                    pad_disp = fe_model.dislocation_dict(i).getDisplacementAtAtom(refpadpos(j,:), fe_model.a, fe_model.material_dict{1,1});
                    %padpos(j,:) = padpos(j,:) + pad_disp;
                    additional_displacements_pad(j,:) = additional_displacements_pad(j,:) + pad_disp;
                end
            end
        end
        
        if length(fe_model.dislocation_dict) > 2
            d1_pos(istep) = fe_model.dislocation_dict(3).coordinates(2);
            if length(fe_model.dislocation_dict) > 3
                d2_pos(istep) = fe_model.dislocation_dict(4).coordinates(2);
            end
        end
    end
    
    if 0>1
        dump = getPosDump(example_name, mstep, refpos, pos, fe_model, fe_analysis, istep, dump);
    end
    
    if trigger == 1
        inv_cont_for = settings.cont_delay; %refresh no_check
        if settings.trigger_case ~= 0 %if settings.trigger_case == 0, this isnt needed because the continuum is always on
            if cont_is_off
                cont_is_off = false;
                cont_on(istep)=istep;
            end
            trigger_hist(istep) = istep;
        end
    else
        inv_cont_for = inv_cont_for -1;
        if settings.trigger_case ~= 0
            if inv_cont_for == 0
                cont_is_off = true;
                cont_off(istep) = istep;
                last_off = istep;
            end
        end
    end
    
    if settings.recover_velocity_center_of_mass
        vel = recoverVelocityCenterOfMass(vel, mxatm);
    end
    
    etot=ekin+epot;
    
    if istep>1
        curr_mean_temp = (curr_mean_temp * (istep-1) + temp)/istep;
    else
        curr_mean_temp = temp;
    end
    
    depth = settings.damping_band_depth;
    minx = min(refpos(:,1));
    miny = min(refpos(:,2));
    maxx = max(refpos(:,1));
    maxy = max(refpos(:,2));
    
    velcount_inner = 0;
    velcount_outer = 0;
    sumv2_inner = 0;
    sumv2_outer = 0;
    
    atoms = [1:mxatm];
    atoms(ismember(atoms,inactive_atoms))=[];
    atoms(ismember(atoms,exclude_atoms))=[];
    for ia=atoms
        if refpos(ia,1) > minx+depth && refpos(ia,1) <= minx+depth+20*a
            velcount_outer = velcount_outer + 1;
            sumv2_outer=sumv2_outer+vel(ia,1)^2+vel(ia,2)^2;
        elseif refpos(ia,1) > minx+depth+20*a
            velcount_inner = velcount_inner + 1;
            sumv2_inner=sumv2_inner+vel(ia,1)^2+vel(ia,2)^2;
        end
    end
    temp_inner = amass/(kb*(2*velcount_inner))*sumv2_inner;
    temp_outer = amass/(kb*(2*velcount_outer))*sumv2_outer;
    temp_inner_container(istep) = temp_inner;
    temp_outer_container(istep) = temp_outer;
    
    tout = temp;
    if strcmp(example_name, 'Tensile_Test') || strcmp(example_name, 'Tensile_Test_rotated')
        if ~mod(istep,100)
            msg_out = sprintf('(%i steps): Sigma_xx = %5.2f Area: %1.5e v1m: %5.2f Strain1: %5.2f Strain2: %5.2f Temp: %5.2f K.\n',istep,dump.mean_stress_xx(istep),stress_area,v1m,strain_hist(istep)*100,strain_hist2(istep)*100, tout);
            fprintf(msg_out)
        end
    elseif strcmp(example_name, 'Example2_rotated') 
        if ~mod(istep,100)
        curr_disp = interface_position_hist(istep) - mean(refpos(interface_atom_list,1));
        target_disp = -495/995*settings.continuum_external_displacement;
        delta_disp = (curr_disp-target_disp)/target_disp*100;
        msg_out = sprintf('(%i steps): Delta_disp = %5.2f per cent, Temp: %5.2f, Temp Interface: %5.2f, Sigma_xx = %5.2f RMSD = %1.5e, (ref = %1.5e)\n',istep,delta_disp,temp,temp_outer,dump.mean_stress_xx(istep),dump.RMSD(istep),max(dump.RMSD_ref));
        fprintf(msg_out)
        end
    end
    
    temp_container(istep) = temp;
    mean_temp_container(istep) = curr_mean_temp;
    temp2_container(istep) = temp^2;
    ekin_container(istep) = ekin;
    epot_container(istep) = epot;
    energy_container(istep) = epot;
    
    if ~mod(istep,iplot)
        if settings.plot_positions
            set(plotatoms,'xData', pos(:,1))
            set(plotatoms,'yData', pos(:,2))
            set(plotpadatoms,'xData', padpos(:,1))
            set(plotpadatoms,'yData', padpos(:,2))
            delete(plotdisloc);
            for i = fe_model.dislocation_dict
                plotdisloc(i.number) = scatter(i.coordinates(1),i.coordinates(2),'o','MarkerEdgeColor',[0 1 0],'MarkerFaceColor',[0 1 0],'SizeData',30);
            end
            title([' Timestep: ', num2str(istep)])
            %drawnow
        end
        set(plotenergy,'xData',1:istep)
        set(plotenergy,'yData',energy_container)
        set(plottemp,'xData',1:istep)
        set(plottemp,'yData',temp_container)
        set(plotmeantemp,'xData',1:istep)
        set(plotmeantemp,'yData',mean_temp_container)
        drawnow
    end
    
    if isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Analysis')
        if ~mod(istep,settings.delta_t_static)
            avgcount = avgcount+1;
            for j = 1:length(interface_atom_list)
                xtsample(avgcount,j) = pos(interface_atom_list(j),1);
                ytsample(avgcount,j) = pos(interface_atom_list(j),2);
            end
            xtmean = mean(xtsample);
            ytmean = mean(ytsample);
            avgcount = 0;
            
             if inv_cont_for > -1
                 if ~mod(istep, irunDD)
                     updateDD = 1;
                 end
                 write_to_nodes = 1;
                 
                 i_entry=0;
                 
                 for i = 1:length(interface_atom_list)
                     for j = 1:fe_model.dimension
                         i_entry=i_entry+1;
                         global_col_index = (interface_node_list(i)-1)*fe_model.dimension + j;
                         boundary_cond_matrix_interface(i_entry,1) = global_col_index;
                         if j==1
                             boundary_cond_matrix_interface(i_entry,2) = (xtmean(i)-refpos(interface_atom_list(i),j));
                         else
                             boundary_cond_matrix_interface(i_entry,2) = (ytmean(i)-refpos(interface_atom_list(i),j));
                         end
                         boundary_cond_matrix_interface(i_entry,3) = vel(interface_atom_list(i),j); %not used in static continuum
                         boundary_cond_matrix_interface(i_entry,4) = force(interface_atom_list(i),j)/amass; %not used in static continuum
                         boundary_cond_matrix_interface(i_entry,5) = interface_node_list(i);
                         boundary_cond_matrix_interface(i_entry,6) = j;
                     end
                 end
                 
                 fe_step = fe_step+1;
                 paddisp = fe_analysis.run(write_to_nodes, boundary_cond_matrix_interface, fe_step, [refpos;refpadpos], [pos; padpos], mxatm, mxpatm, updateDD, istep, settings);
                 write_to_nodes = 0;
                 updateDD = 0;
                 
                 padpos = refpadpos + additional_displacements_pad + paddisp;
                 
                 avg_counter = 0;
             end
             
             
        else
            avgcount = avgcount+1;
            for j = 1:length(interface_atom_list)
                xtsample(avgcount,j) = pos(interface_atom_list(j),1);
                ytsample(avgcount,j) = pos(interface_atom_list(j),2);
            end
        end
        
    elseif isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Dynamic_Analysis')
        if inv_cont_for > -1
            if ~mod(istep, iwriteVTK)
                write_to_nodes = 1;
            end
            
            if ~mod(istep, irunDD)
                updateDD = 1;
            end
            
            i_entry=0;
            
            for i = 1:length(interface_atom_list)
                for j = 1:fe_model.dimension
                    i_entry=i_entry+1;
                    global_col_index = (interface_node_list(i)-1)*fe_model.dimension + j;
                    boundary_cond_matrix_interface(i_entry,1) = global_col_index;
                    boundary_cond_matrix_interface(i_entry,2) = (pos(interface_atom_list(i),j)-refpos(interface_atom_list(i),j));
                    boundary_cond_matrix_interface(i_entry,3) = vel(interface_atom_list(i),j);
                    boundary_cond_matrix_interface(i_entry,4) = force(interface_atom_list(i),j)/amass;
                    boundary_cond_matrix_interface(i_entry,5) = interface_node_list(i);
                    boundary_cond_matrix_interface(i_entry,6) = j;
                end
            end
            
            fe_step = fe_step+1;
            paddisp = fe_analysis.run(write_to_nodes, boundary_cond_matrix_interface, fe_step, [refpos;refpadpos], [pos; padpos], mxatm, mxpatm, updateDD, istep, settings);
            write_to_nodes = 0;
            updateDD = 0;
            
            padpos = refpadpos + additional_displacements_pad + paddisp;
            
            avg_counter = 0;
        end
        
    elseif isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Hybrid_Analysis')
        write_to_nodes = 0;
        
        %DYNAMIC PART
        if inv_cont_for > -1 || istep >= settings.continuum_external_change(1) && istep <= settings.continuum_external_change(2)
            if ~mod(istep, irunDD)
                updateDD = 1;
            end
            i_entry=0;
            for i = 1:length(interface_atom_list)
                for j = 1:fe_model.dimension
                    i_entry=i_entry+1;
                    global_col_index = (interface_node_list(i)-1)*fe_model.dimension + j;
                    boundary_cond_matrix_interface(i_entry,1) = global_col_index;
                    boundary_cond_matrix_interface(i_entry,2) = (pos(interface_atom_list(i),j)-refpos(interface_atom_list(i),j));
                    boundary_cond_matrix_interface(i_entry,3) = vel(interface_atom_list(i),j);
                    boundary_cond_matrix_interface(i_entry,4) = force(interface_atom_list(i),j)/amass;
                    boundary_cond_matrix_interface(i_entry,5) = interface_node_list(i);
                    boundary_cond_matrix_interface(i_entry,6) = j;
                end
            end
            fe_step = fe_step+1;
            paddisp_dyn = fe_analysis.dynamic_analysis.run(write_to_nodes, boundary_cond_matrix_interface, fe_step, [refpos;refpadpos], [pos; padpos], mxatm, mxpatm, updateDD, istep, settings);
            updateDD = 0;
        end
        
        %STATIC PART
        if ~mod(istep,settings.delta_t_static) && istep >= settings.continuum_external_change(1) && istep <= settings.continuum_external_change(2)
            i_entry=0;
            %set interface BC to zero
            for i = 1:length(interface_atom_list)
                for j = 1:fe_model.dimension
                    i_entry=i_entry+1;
                    global_col_index = (interface_node_list(i)-1)*fe_model.dimension + j;
                    boundary_cond_matrix_interface(i_entry,1) = global_col_index;
                    boundary_cond_matrix_interface(i_entry,2) = 0;
                    boundary_cond_matrix_interface(i_entry,3) = 0;
                    boundary_cond_matrix_interface(i_entry,4) = 0;
                    boundary_cond_matrix_interface(i_entry,5) = interface_node_list(i);
                    boundary_cond_matrix_interface(i_entry,6) = j;
                end
            end
            fe_step = fe_step+1;
            paddisp_stat = fe_analysis.static_analysis.run(write_to_nodes, boundary_cond_matrix_interface, fe_step, [refpos;refpadpos], [pos; padpos], mxatm, mxpatm, updateDD, istep, settings);
            
        end
        
        if ~mod(istep, iwriteVTK)
            write_to_nodes = 1;
            %compute u_dd
            if write_to_nodes
                
                %TODO: DD part not correct
                total_displacement = fe_analysis.dynamic_analysis.u + fe_analysis.static_analysis.u;                
                if ~isempty(fe_analysis.dynamic_analysis.model.dislocation_dict)
                    for node = fe_analysis.dynamic_analysis.model.node_dict
                        for coord_index = 1:fe_analysis.dynamic_analysis.model.dimension
                            global_index = (node.number-1) * fe_analysis.dynamic_analysis.model.dimension + coord_index;
                            fe_analysis.dynamic_analysis.u_dd(global_index) = self.getDDDisplacementatNode(node, coord_index);
                        end
                    end
                end
                fe_analysis.dynamic_analysis.updateDOF(total_displacement + fe_analysis.dynamic_analysis.u_dd, total_displacement)
            end
        end
        
        padpos = refpadpos + additional_displacements_pad;
        if exist('paddisp_dyn', 'var')
            padpos = padpos + paddisp_dyn;
        end
        if exist('paddisp_stat', 'var')
            padpos = padpos + paddisp_stat;
        end
    end
    
    if ~mod(istep, inl) %Update Neighborlist
        [ nl, numneigh] = neighborlist([pos; padpos], mxatm, mxpatm, rcut, pbc, xbox, ybox, a, example_name, settings);
    end
    
    if ~mod(istep, iwriteVTK)
        %Write ouput to vtk file
        fe_analysis.output_handler.write(fe_model, vtks, [refpos; refpadpos], [pos; padpos], mxatm, mxpatm, 0, energies, vel, fe_analysis.u1);
        vtks = vtks + 1;
    end
    
    if istep == mstep_save_at
        disp('Saving current state of simulation to file ...')
        mkdir(['examples/',fe_analysis.example,'/restart_sim'])
        save(['examples/',fe_analysis.example,'/restart_sim/workspace_',cont_type,'_',num2str(tset),'K.mat'],'-v7.3')
        disp('... Done saving')
    end
    
end

if settings.trigger_case ~= 0
    if length(nonzeros(cont_off)) ~= length(nonzeros(cont_on))
        cont_off(mstep) = mstep;
    end
end

disp('fin')


if settings.write_reflection_file
reflected3 = max(max_disp_x_hist(settings.int_back(1):settings.int_back(2)))/max(max_disp_x_hist(settings.int_forth(1):settings.int_forth(2)))*100;

fileID = fopen([fe_analysis.output_handler.folder_name,'/reflection.txt'],'w');
fprintf(fileID,'reflected3 %f \n',reflected3);
fclose(fileID);

end

if settings.save_workspace_after_finishing
    if isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Analysis') || isa(fe_analysis,'nsAnalyzer.nsAnalysis.Linear_Hybrid_Analysis')
        disp(['Saving workspace to ',fe_analysis.output_handler.file_name,'/',sim_name,'_',num2str(settings.delta_t_static),'.mat ...'])
        mkdir([fe_analysis.output_handler.file_name])
        delete([fe_analysis.output_handler.file_name,'/*.mat'])
        save([fe_analysis.output_handler.file_name,'/',sim_name,'_',num2str(settings.delta_t_static),'.mat'],'-v7.3')
    else
        disp(['Saving workspace to ',fe_analysis.output_handler.file_name,'/',sim_name,'.mat ...'])
        mkdir([fe_analysis.output_handler.file_name])
        delete([fe_analysis.output_handler.file_name,'/*.mat'])
        save([fe_analysis.output_handler.file_name,'/',sim_name,'.mat'],'-v7.3')
    end
    disp('... Done saving')
end
end

