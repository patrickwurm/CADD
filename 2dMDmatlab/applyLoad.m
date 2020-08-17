function pos = applyLoad(inactive_atoms, atom_load, pos, refpos, example_name, istep, settings)


if strcmp(example_name, 'Example2') || strcmp(example_name, 'Example2_rotated') 
    
    pulse_step_1 = atom_load.pulse_step_1;
    period_1 = atom_load.period_1;
    amplitude_1 = atom_load.amplitude_1;
    pulse_step_2 = atom_load.pulse_step_2;
    period_2 = atom_load.period_2;
    amplitude_2 = atom_load.amplitude_2;
    pulse_step_3 = atom_load.pulse_step_3;
    period_3 = atom_load.period_3;
    amplitude_3 = atom_load.amplitude_3;
    pulse_step_4 = atom_load.pulse_step_4;
    period_4 = atom_load.period_4;
    amplitude_4 = atom_load.amplitude_4;
    pulse_step_5 = atom_load.pulse_step_5;
    period_5 = atom_load.period_5;
    amplitude_5 = atom_load.amplitude_5;
    pulse_step_6 = atom_load.pulse_step_6;
    period_6_1 = atom_load.period_6_1;
    period_6_2 = atom_load.period_6_2;
    period_6_3 = atom_load.period_6_3;
    amplitude_6_1 = atom_load.amplitude_6_1;
    amplitude_6_2 = atom_load.amplitude_6_2;
    amplitude_6_3 = atom_load.amplitude_6_3;
    pulse_step_7 = atom_load.pulse_step_7;
    period_7 = atom_load.period_7;
    amplitude_7 = atom_load.amplitude_7;
    
%     for i = 1:length(inactive_atoms)
%         if (istep>=pulse_step_1 && istep<pulse_step_1+period_1)
%             pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_1*sin((istep-pulse_step_1)/period_1*2*pi/2);
%         end
%         %             if (istep>=pulse_step_2 && istep<pulse_step_2+period_2)
%         %                 pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_2*sin((istep-pulse_step_2)/period_2*2*pi/2);
%         %             end
%         if (istep>=pulse_step_3 && istep<pulse_step_3+period_3)
%             pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_3*sin((istep-pulse_step_3)/period_3*2*pi/2);
%         end
%         if (istep>=pulse_step_4 && istep<pulse_step_4+period_4)
%             pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_4*sin((istep-pulse_step_4)/period_4*2*pi/2);
%         end
%         if (istep>=pulse_step_5 && istep<pulse_step_5+period_5)
%             pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_5*sin((istep-pulse_step_5)/period_5*2*pi/2);
%         end
%         if (istep>=pulse_step_6 && istep<pulse_step_6+period_6_1)
%             pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_6_1*sin((istep-pulse_step_6)/period_6_1*2*pi/2)-amplitude_6_2*sin((istep-pulse_step_6)/period_6_2*2*pi/2)-amplitude_6_3*sin((istep-pulse_step_6)/period_6_3*2*pi/2);
%         end
%         if (istep>=pulse_step_7 && istep<pulse_step_7+period_7)
%             pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_7*sin((istep-pulse_step_7)/period_7*pi/2);
%         end
%     end

if 0>1
    for i = 1:length(inactive_atoms)
        if (istep>=pulse_step_1 && istep<pulse_step_1+period_1)
            pos(inactive_atoms(i),1) = refpos(inactive_atoms(i),1) - amplitude_1*sin((istep-pulse_step_1)/period_1*2*pi/2);
        end
    end
end

if isfield(settings.additional_input,'pulse')
    if settings.additional_input.pulse > 0
        n_init = 1000;
        alpha=0.05;
        omega = 1/period_1*2*pi;
        preruntime = 3*sqrt(1/(alpha*omega^2));
        time = istep-n_init-preruntime;
        psi = real( 1/2*amplitude_1 * exp(-alpha*(omega*time)^2) * ( exp(1j*omega*time) + exp(-1j*omega*time) ));
        pos(inactive_atoms,1) = refpos(inactive_atoms,1) - psi;
    end
end
    
elseif strcmp(example_name, 'Dislocation')
    
    if istep > indent_start && istep <= indent_start + indent_duration
        pos(inactive_atoms,2)=refpos(inactive_atoms,2)-indent_depth/indent_duration*(istep-indent_start);
        pos(inactive_atoms,1)=refpos(inactive_atoms,1);
    elseif istep > indent_start + indent_duration
        pos(inactive_atoms,2) = refpos(inactive_atoms,2)-indent_depth;
        pos(inactive_atoms,1)=refpos(inactive_atoms,1);
    end
    
elseif strcmp(example_name, 'Radial_Pulse')
    
    midx = atom_load.midx;
    midy = atom_load.midy;
    Amp = atom_load.Amp;
    sig = atom_load.sig;
    uc = atom_load.uc;
    a = atom_load.a;
    indent_start = atom_load.indent_start;
    
    if istep == indent_start
        for k = 1:length(inactive_atoms)
            r = sqrt( (refpos(inactive_atoms(k),1)-midx)^2 + (refpos(inactive_atoms(k),2)-midy)^2 );
            angle = atan( (refpos(inactive_atoms(k),2)-midy) / (refpos(inactive_atoms(k),1)-midx) );
            rdisp = Amp/(Amp-uc)*(Amp*exp(-(r/sig)^2)-uc);
            
            if refpos(inactive_atoms(k),1)-midx < 0 - 0.01*a
                pos(inactive_atoms(k),1) = pos(inactive_atoms(k),1) - rdisp * cos(angle);
                pos(inactive_atoms(k),2) = pos(inactive_atoms(k),2) - rdisp * sin(angle);
            else
                pos(inactive_atoms(k),1) = pos(inactive_atoms(k),1) + rdisp * cos(angle);
                pos(inactive_atoms(k),2) = pos(inactive_atoms(k),2) + rdisp * sin(angle);
            end
        end
    end
elseif strcmp(example_name, 'Tensile_Test') || strcmp(example_name, 'Tensile_Test_rotated')
    %do nothing
    
else
    error('Example not found.')
end
    
end

