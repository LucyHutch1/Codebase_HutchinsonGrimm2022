%% Oliver Grimm and Lucy Hutchinson 2021
%
function test_model_predictions(filename_testdata, filename_pat_info, filename_pop_par, directory,task_name)
%% Test the model predictions
% Load the script that summarises the data for the holdout pats
load(filename_testdata);
load(filename_pat_info);

pats_chosen = data_summary.PatID;
endtimes_curr = data_summary.endTimes;
tilenumbers = data_summary.Ntilespre;
tile_limit =50;
n_stoch = 5;


% Now load the pop parameters
load(filename_pop_par);

%% 1. Run the model for all of the holdout patients using the derived
% population parameters saving every timepoint. Also run for the sets of
% control parameters for comparison

% RUN POP PARS
par_names = {'IMpprol','IMpkill','IMpdeath','IMrwalk'};
parsets = pop_param_vals;
n_pars = length(parsets);
% pass in endtimes+10
sim_endtimes = endtimes_curr+10;
 
 
% Make a new directory to save the results
mkdir(directory);
save_path = [directory,'/'];
savetimecourse = true;
saveselectedgridstack = nan(length(pats_chosen),4);
for p = 1:length(pats_chosen)
saveselectedgridstack(p,:) = [pats_chosen(p),randi(tilenumbers(p)), randi(n_pars), randi(n_stoch)];
end

% Run the simulations with pop params
runsims_parsets(parsets, par_names, pats_chosen, endtimes_curr, sim_endtimes,...
   tilenumbers, tile_limit,  n_stoch, save_path, task_name,savetimecourse,saveselectedgridstack)%


% RUN CONTROL PARS
savetimecourse = true;
saveselectedgridstack = [];
 for ii =1:5 %hardcoded 5- consider passing in
parsets = control_param_vals(:,:,ii);
task_name = ['control',num2str(ii)];
runsims_parsets(parsets, par_names, pats_chosen, endtimes_curr, sim_endtimes,...
    tilenumbers, tile_limit,  n_stoch, save_path, task_name,savetimecourse,saveselectedgridstack);
 end

