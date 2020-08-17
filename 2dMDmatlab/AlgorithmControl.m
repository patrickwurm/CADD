function [inv_cont_for, trigger, run_in, pos_matrix, band_atom_ass_atoms_position, filter_output] = AlgorithmControl(trigger_case, pos_matrix, kernel_matrix, band_atoms, band_atom_ass_atoms, band_atom_ass_atoms_position, run_in, mstep_run_in, example, load_fluc, restart_sim, D, istep, inv_cont_for, pos)

filter_output = [];

if istep == mstep_run_in+1
    if load_fluc == true
        load(['examples/',example,'/run_in_results/recorded_eq_fluc.mat'],'run_in')
        inv_cont_for = 0;
        trigger = 0;
        max_distance_old = run_in.max_distance;
        min_distance_old = run_in.min_distance;
        run_in.max_distance = run_in.max_distance; % + (max_distance_old - min_distance_old)*0.1;
        run_in.min_distance = run_in.min_distance; % - (max_distance_old - min_distance_old)*0.1;
        clear max_distance_old min_distance_old
    end
end
if istep == mstep_run_in && load_fluc == false %&& restart_sim ~= 1
    run_in.max_distance = max(run_in.rel_distance,[],3);
    run_in.min_distance = min(run_in.rel_distance,[],3);
    promptMessage = sprintf('Saving recorded equilibrium fluctuations. Are you sure you want to overwrite existing data?');
    button = questdlg(promptMessage, 'Confirm overwrite', 'Yes', 'No', 'Cancel', 'Yes');
    if strcmp(button,'Yes')
    disp(['Saving recorded equilibrium fluctuations at time step ',num2str(mstep_run_in), ' (mstep_run_in) to file examples/',example,'/run_in_results/recorded_eq_fluc.mat',])
    save(['examples/',example,'/run_in_results/recorded_eq_fluc.mat'],'run_in')
    inv_cont_for = 0;
    trigger = 0;
    elseif strcmp(button,'No')
        disp(['I did not write the recorded equilibrium fluctuations to file.'])
    elseif strcmp(button,'Cancel')
        error('Stopping execution now.')
    end
end

if istep > mstep_run_in
    if ~mod(istep,D)
        i=0;
        for j = band_atom_ass_atoms
            i = i+1;
            band_atom_ass_atoms_position(i,:) = pos(j,:);
        end
        
        pos_matrix = cat(3, [pos(band_atoms,1)-pos(band_atom_ass_atoms,1), pos(band_atoms,2)-pos(band_atom_ass_atoms,2)], pos_matrix(:,:,1:end-1));
        
        [trigger, filter_output] = Cont_trigger_DSP(trigger_case, 1, pos_matrix, kernel_matrix, band_atom_ass_atoms_position, run_in);
        
        %if trigger == 1
        %    disp(['Activation at step: ', num2str(istep)])
        %end
    else
        trigger=0;
    end
elseif istep <= mstep_run_in
    if ~mod(istep,D)
        i=0;
        for j = band_atom_ass_atoms
            i = i+1;
            band_atom_ass_atoms_position(i,:) = pos(j,:);
        end
        pos_matrix = cat(3, [pos(band_atoms,1)-pos(band_atom_ass_atoms,1), pos(band_atoms,2)-pos(band_atom_ass_atoms,2)], pos_matrix(:,:,1:end-1));
        
        [trigger, filter_output] = Cont_trigger_DSP(trigger_case, 0, pos_matrix, kernel_matrix, band_atom_ass_atoms_position, run_in);
        run_in.rel_distance(:,:,istep/D) = filter_output;
    else
        trigger=1;
    end
end


end

