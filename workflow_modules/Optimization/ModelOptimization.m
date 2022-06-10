%% Oliver Grimm and Lucy Hutchinson 2021
%
function ModelOptimization(n_pars, par_names,filenameout_pars, pats_chosen, datasummary_to_load, directory, save_path, task_name)
%% Model Optimization
% This function generates parameter sets and simulates all patients, all
% tiles in the training set with all parameter configurations in order to investigate
% which parameter configurations produce the best reproduction of the
% observed on-treatment samples

%% 1. Generate samples

% Ranges:
% IMpprol 0.01 - 0.4 (log)
% IMpkill 0.001 - 0.5 (log)
% IMpdeath 0.001 - 0.2 (log)
% IMrwalk 0.7-1

% Use the following line to generate new parameter samples
pars_all = [10.^(-2+rand([n_pars,1])*1.6021), 10.^(-3+rand([n_pars,1])*2.699),10.^(-3+rand([n_pars,1])*2.3010), 0.7+rand([n_pars,1])*0.3];

%% 2. Run simulations
% Load patient data
load(datasummary_to_load);

if isempty(pats_chosen)
pats_chosen = data_summary.PatID;% Run the optimization for all patients in the training set
end

tile_number_patschosen=[];endtimes_patschosen=[];

for ii = 1:length(pats_chosen)
tile_number_patschosen(ii) = data_summary.Ntilespre(data_summary.PatID==pats_chosen(ii)); % Get the tile numbers so it is easily accessible
endtimes_patschosen(ii) = data_summary.endTimes(data_summary.PatID==pats_chosen(ii)); % Get the endtimes as well
end
tile_limit=30;

n_stoch = 5; % stochastic repeats for each tile, each parameter set
tile_inds_all = {}; % Save the actual tiles simulated for each patient

mkdir(directory);
% Save pars in new directory
save(filenameout_pars,'pars_all','par_names'); 

savetimecourse = false;
saveselectedgridstack = [];

runsims_parsets(pars_all, par_names, pats_chosen, endtimes_patschosen, endtimes_patschosen,...
    tile_number_patschosen, tile_limit,  n_stoch, save_path, task_name,savetimecourse,saveselectedgridstack)%
