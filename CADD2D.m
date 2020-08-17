function CADD2D(example, cont_type, additional_input)

close all
clearvars -except example cont_type additional_input
restoredefaultpath

addpath([cd,'/soofeaM_combined_explicit_dynamic/src/']);
addpath([cd,'/examples/',example]);
addpath([cd,'/2dMDmatlab/']);

fprintf('I am creating the FE Model ...\n')
output_file_name = example;

if strcmp(cont_type,'static') || strcmp(cont_type,'hybrid')
    if strcmp(additional_input.example_id_1,'Reflection')
        if additional_input.pulse > 0
            inp_vec_0 = [34];
        elseif additional_input.pulse == 0
            inp_vec_0 = 35;
        end
    else
        inp_vec_0 = 31;
    end
else
    inp_vec_0 = 10^6; %arbitrary
end

for i = 1:length(inp_vec_0)
    
    if isfield(additional_input,'trigger_case')
        settings.trigger_case = additional_input.trigger_case; %0: CADD, 7: CADD+A, 100: rigid continuum
    else
        settings.trigger_case = 0; %0: CADD, 7: CADD+A, 100: rigid continuum
    end
    settings.cont_delay = 2500; %Continuum delay, how many time steps must the continuum be active after activation 
    settings.bandatom_position = 50; %distance of band of atoms from the interface
    settings.delta_t_static = inp_vec_0(i); %continuum time step used in the quasi-static continuum
    
    sim_name = cont_type;
    
    if strcmp(additional_input.example_id_1,'Reflection')
        pulse_number = additional_input.pulse;
        prepath = ['examples/',example,'/Reflection/Puls',num2str(pulse_number)];
    elseif strcmp(additional_input.example_id_1,'Tensile_Test') 
        prepath = ['examples/',example,'/Tensile_Test/',num2str(additional_input.temperature),'K/',additional_input.example_id_2,'/Seed_',num2str(additional_input.rng_seed)];
    elseif strcmp(additional_input.example_id_1,'Equilibrium') || strcmp(additional_input.example_id_1,'Equilibrium_small')
        prepath = ['examples/',example,'/',additional_input.example_id_1,'/',num2str(additional_input.temperature),'K/Speed',num2str(additional_input.speed),'/',additional_input.example_id_2];
    end

    if strcmp(cont_type,'static')
        results_path = [prepath,'/results_stat_',num2str(inp_vec_0(i)),'/results'];
        outfile_path = [results_path,'/output'];
    elseif strcmp(cont_type,'hybrid')
        results_path = [prepath,'/results_hyb_',num2str(inp_vec_0(i)),'/results'];
        outfile_path = [results_path,'/output'];
    elseif strcmp(cont_type,'dynamic')
        results_path = [prepath,'/results_dyn/results'];
        outfile_path = [results_path,'/output'];
    end
    mkdir(outfile_path)
    settings.additional_input = additional_input;
    [model, analysis] = feval(example,outfile_path, results_path, cont_type, settings);
        
    TwoDMDmatlab(model, analysis, settings, sim_name)
    
end

end
