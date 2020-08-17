function force_tables = initializeRainbowTables(a, rcut, atomic_potential)

if atomic_potential ~= 4
    error('Rainbow tables currently only implemented for EAM potential')
end

convA_to_meters = 10^-10;
conveV_to_J = 1.602176565*10^-19;

%ZHOU 2001
% re = 2.886166*convA_to_meters;
% fe = 1.392302*conveV_to_J/convA_to_meters;
% rhoe = 20.226537*conveV_to_J/convA_to_meters;
% rhos = rhoe;
% alpha = 6.942419;
% beta = 3.702623;
% A = 0.251519*conveV_to_J;
% B = 0.313394*conveV_to_J;
% kappa = 0.395132;
% lambda = 0.790264;
% Fn = [-2.806783; -0.276173; 0.893409; -1.637201]*conveV_to_J;
% F = [-2.83; 0; 0.929508; -0.682320]*conveV_to_J;
% eta = 0.779208;
% Fe = -2.829437*conveV_to_J;

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

dr = 2001;
dx = rcut/(dr-1);

force_tables.X = [0:dx:rcut+dx]';
force_tables.rhoi = fe*exp(-beta*(force_tables.X/re-1))./ (1+(force_tables.X/re-lambda).^20);

force_tables.pi1 = fe*exp(-beta*(force_tables.X/re-1))*(-beta/re).*(1+(force_tables.X/re-lambda).^20).^-1 ...
        -fe*exp(-beta*(force_tables.X/re-1)).*(1+(force_tables.X/re-lambda).^20).^-2*20.*(force_tables.X/re-lambda).^19.*1/re;

force_tables.phi = A*exp(-alpha*(force_tables.X/re-1))./(1+(force_tables.X/re-kappa).^20)...
        -B*exp(-beta*(force_tables.X/re-1))./(1+(force_tables.X/re-lambda).^20);
    
force_tables.phi1 = A*exp(-alpha*(force_tables.X/re-1))*(-alpha/re).*(1+(force_tables.X/re-kappa).^20).^(-1) ...
        - A*exp(-alpha*(force_tables.X/re-1)).*(1+(force_tables.X/re-kappa).^20).^(-2)*20.*(force_tables.X/re-kappa).^19*1/re ...
        - B*exp(-beta*(force_tables.X/re-1))*(-beta/re).*(1+(force_tables.X/re-lambda).^20).^(-1) ...
        + B*exp(-beta*(force_tables.X/re-1)).*(1+(force_tables.X/re-lambda).^20).^(-2)*20.*(force_tables.X/re-lambda).^19*1/re;                

force_tables.use_them = 1;

%truncate, so that entries at rcut = 0;
force_tables.rhoi = force_tables.rhoi - force_tables.rhoi(end-1);
force_tables.pi1 = force_tables.pi1 - force_tables.pi1(end-1);
force_tables.phi = force_tables.phi - force_tables.phi(end-1);
force_tables.phi1 = force_tables.phi1 - force_tables.phi1(end-1);

%set last entry also to zero, this entry is needed if 
%rrr is slightly larger than rcut due to round-off error
force_tables.rhoi(end) = 0; 
force_tables.pi1(end) = 0; 
force_tables.phi(end) = 0; 
force_tables.phi1(end) = 0; 
end

