function atom_load = setupLoad(example_name, a, axmin, axmax, aymin, aymax, rcut,settings)

if strcmp(example_name, 'Example2') || strcmp(example_name, 'Example2_rotated')
    
    atom_load.pulse_step_1=1000;
    atom_load.period_1 = 10^6;
    if isfield(settings.additional_input,'pulse')
        if settings.additional_input.pulse == 1
            atom_load.period_1 = 2500;
        elseif settings.additional_input.pulse == 2
            atom_load.period_1 = 5000;
        elseif settings.additional_input.pulse == 3
            atom_load.period_1 = 7500;
        elseif settings.additional_input.pulse == 0
            atom_load.period_1 = 666; %arbitrary, not used
        end
    end
    %atom_load.period_1 = 1500; %Solitary 1
    atom_load.amplitude_1 = 0.25*a; %Solitary 2
    %atom_load.amplitude_1 = 0.025*a; %Solitary 1 
    atom_load.pulse_step_2=200000;
    atom_load.period_2 = 1000;
    atom_load.amplitude_2 = 2*a;
    atom_load.pulse_step_3=200000;
    atom_load.period_3 = 10000;
    atom_load.amplitude_3 = 2*a;
    atom_load.pulse_step_4=400000;
    atom_load.period_4 = 50000;
    atom_load.amplitude_4 = 2*a;
    atom_load.pulse_step_5=500000;
    atom_load.period_5 = 50000;
    atom_load.amplitude_5 = 4*a;
    atom_load.pulse_step_6=600000;
    atom_load.period_6_1 = 50000;
    atom_load.period_6_2 = 10000;
    atom_load.period_6_3 = 5000;
    atom_load.amplitude_6_1 = 2*a;
    atom_load.amplitude_6_2 = 1*a;
    atom_load.amplitude_6_3 = 0.5*a;
    atom_load.pulse_step_7=900000;
    atom_load.period_7 = 50000;
    atom_load.amplitude_7 = 4*a;
    
elseif strcmp(example_name, 'Dislocation')
    
    atom_load.a = a;
    atom_load.indent_start = 1000;
    atom_load.indent_depth = 3*a;
    atom_load.indent_duration = 10000;
    
elseif strcmp(example_name, 'Radial_Pulse')
        
    atom_load.midx = (axmax+axmin)/2;
    atom_load.midy = (aymax+aymin)/2;
    atom_load.Amp = 2*10^-10;
    atom_load.sig = 3*a;
    atom_load.uc = atom_load.Amp*exp(-(rcut/atom_load.sig)^2);
    atom_load.a = a;
    atom_load.indent_start = 1000;
    atom_load.indent_depth = 3*a;
    atom_load.indent_duration = 10000;

        
elseif strcmp(example_name, 'Tensile_Test') || strcmp(example_name, 'Tensile_Test_rotated')
        
    atom_load = [];
    
else
    error('Example not implemented.')
end

end

