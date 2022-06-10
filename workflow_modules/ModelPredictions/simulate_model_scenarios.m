%% Oliver Grimm and Lucy Hutchinson 2021
%
function simulate_model_scenarios(filename_test_data_summary, filename_patinfo, filename_pop_par, directory,task_name)

% Load the script that summarises the data for the holdout pats
load(filename_test_data_summary);
load(filename_patinfo);

pats_chosen = data_summary.PatID;
tilenumbers = data_summary.Ntilespre;
tile_limit =50;
n_stoch = 5;


% Now load the pop parameters
load(filename_pop_par);

%% 1. Run the model for all of the holdout patients using the derived
% population parameters and two variations to represent different treatments saving every timepoint. 

% Scenario 1: pop pars. 
% Scenario 2: pop pars with higher IMpprol *1.5
% Scenario 3. pop pars with higher IMinflRate *8
par_names_all{1} = {'IMpprol','IMpkill','IMpdeath','IMrwalk'};
par_names_all{2} = {'IMpprol','IMpkill','IMpdeath','IMrwalk'};
par_names_all{3} = {'IMpprol','IMpkill','IMpdeath','IMrwalk','IMinflRate'};

parsets{1} = pop_param_vals;
parsets{2} = pop_param_vals.*[1.5,1,1,1];
parsets{3} = [pop_param_vals,8*ones(length(pop_param_vals),1)];

n_pars = length(parsets);

% Simulate to day 50 for all pats
sim_endtimes = repmat(50,size(pats_chosen));
 
 
% Make a new directory to save the results
mkdir(directory);
save_path = [directory,'/'];
savetimecourse = true;
saveselectedgridstack = nan(length(pats_chosen),4);
for p = 1:length(pats_chosen)
saveselectedgridstack(p,:) = [pats_chosen(p), randi(n_stoch), 20, 50];% SAVE TIMEPOINTS 20 and 50
end


% Run 3 sets of simulations
for set_ind=1:3
    pars_all = parsets{set_ind};
    par_names = par_names_all{set_ind};
    task_name_in = [task_name,num2str(set_ind)];
% Run the simulations with pop params
runsims_parsets_savemultitile(pars_all, par_names, pats_chosen, endtimes, sim_endtimes,...
    tilenumbers, tile_limit,  n_stoch, save_path, task_name_in,savetimecourse,saveselectedgridstack)%
end