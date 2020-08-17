% adding paths
addpath src
addpath (['examples/',example]);

% Creating the model
fprintf('Creating FE Model ...\n')
[model, analysis] = feval(example);

% running the analysis - whether linear or nonlinear is specified in the
% example file.
%tic
%fprintf('\nStarting Analysis ...\n')
%analysis.run();
%disp('... done')
%toc
