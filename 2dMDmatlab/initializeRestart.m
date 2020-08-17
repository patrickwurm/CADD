function [istep, mstep, iwriteVTK, mstep_save_at, mstep_run_in, restart_sim, i_percent, timer_1, inv_cont_for, load_fluc, settings] = initializeRestart(fe_model, fe_analysis, restart_settings)

fe_analysis.setModel(fe_model);
disp('... Done loading')
disp('-> Resetting timestep')
istep=restart_settings.istep;
mstep=restart_settings.mstep;
iwriteVTK = restart_settings.iwriteVTK;
mstep_save_at = restart_settings.mstep_save_at;
mstep_run_in = restart_settings.mstep_run_in;
restart_sim = restart_settings.restart_sim;
i_percent = restart_settings.i_percent;
timer_1 = tic;
inv_cont_for = restart_settings.inv_cont_for;
load_fluc = restart_settings.load_fluc;
settings = restart_settings.settings;
end

