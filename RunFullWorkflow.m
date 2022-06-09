%% Oliver Grimm and Lucy Hutchinson 2021
%
%% Master script : Biopsy scheduling workflow
% This workflow uses the data provided in the zip folder 

%% 0. Starter script to add all files to the path
RUNFIRST;

%% 1. Prepare a summary table of the processed data
% This script will save two .mat files- one for the training set and one
% for the test set. It also will calculate the RDFs for pre- and post-
% tiles and save them.
%prep_data_summary

%% Parameters applicable for all workflow steps
n_stoch = 5; % stochastic repeats
n_pars =1000;% for Optimization
obs_range_tol = 0.2; % Width of bounds for SAM
frac_within_tol = 0.8; % Fraction of points that must lie within the bounds for SAM acceptance
thresh_acceptSAM = 0.7; % Threshold SAM value for acceptance of a parameter set
thresh_acceptVarSAM = 0.3; % Threshold VarSAM value for acceptance of a parameter set
thresh_votesforaccept = 0.6; % Proportion of patients that must have an acceptable SAM/VarSAM to accept a parameter set
par_names = {'IMpprol','IMpkill','IMpdeath','IMrwalk'}; %
filename_parsets = 'workflow_modules/Optimization/results/parameter_samples_optimization1000.mat';
filename_datasummary = 'final_tiles/PatientData/pat_summary_training.mat';
filename_test_data_summary ='final_tiles/PatientData/pat_summary_test.mat';
filename_patinfo = 'final_tiles/patinfo.mat';

data_summ_train = load('final_tiles/PatientData/pat_summary_training.mat');
data_summ_test = load('final_tiles/PatientData/pat_summary_test.mat');

pats_chosen = data_summ_train.data_summary.PatID;
pats_test = data_summ_test.data_summary.PatID;

%% 2. Sensitivity Analysis
% Runs the full sensitivity analysis and saves a MAT file with
% the results
pats_chosen_SA = [10,14,24,56,33,12]; % In order of increasing number of CD8 cells: ~0,50,100,200,400,800
tiles_chosen_SA = [7,3,6,19,6,5]; % particular tiles selected
directory = ['workflow_modules/SensitivityAnalysis/Results'];
filenameout_SA = ['workflow_modules/SensitivityAnalysis/Results/SAresults_all'];
%Sensitivity_analysis(n_stoch,pats_chosen_SA, tiles_chosen_SA ,directory,filenameout_SA);

% Process and plot the results for the sensitivity analysis
 filename_in =filenameout_SA;
 directory = ['workflow_modules/SensitivityAnalysis/Results/Figures'];
 plotSensitivityAnalysis(filename_in,directory,pats_chosen_SA, tiles_chosen_SA);
% 
% %% 3. Model Optimization: Simulate a large number of parameter sets for all
% % patients in the training set (sampling 30 tiles for those with >30) and
% % saving the relevant results
% directory = 'workflow_modules/Optimization/results';
% save_path = 'workflow_modules/Optimization/results';
% task_name = ('Optimization_1000parsets');
% 
%  ModelOptimization(n_pars,par_names,filename_parsets,pats_chosen, filename_datasummary, directory, save_path, task_name);
% 
% % Analyse optimization results (Get SAM and VarSAM and calculate RDFs)
% results_path = 'workflow_modules/Optimization/results';
% save_path = 'workflow_modules/Optimization/results/processed_optimization_results__';
% analyse_ModelOptimization(obs_range_tol, frac_within_tol, filename_datasummary,pats_chosen, filename_parsets, results_path,task_name, save_path); %
% 
% % Plot the optimization results. In this script we generate heatmaps and diagnostic plots
% opt_results_processed = save_path;
% save_path='workflow_modules/Optimization/results/Figure_'; 
% plot_ModelOptimization(thresh_acceptSAM, thresh_acceptVarSAM, thresh_votesforaccept, obs_range_tol, opt_results_processed,filename_parsets, save_path);
% 
% %% 4. Model testing. Extract pop pars using the SAM and test these pars on unseen patient data
% % Processing of model optimization results to obtain the population
% % parameters. 
% getPopPars(thresh_acceptSAM, thresh_acceptVarSAM,thresh_votesforaccept,...
%       opt_results_processed, filename_parsets,filename_datasummary,results_path, task_name);
% 
% % % Test the model predictions using the holdout patients and population
% % % parameters
%  filename_pop_par = [results_path,'/pop_params'];
%  directory = ['workflow_modules/ModelValidation/Results'];
%  task_name = 'test';
%  test_model_predictions(filename_test_data_summary,filename_patinfo, filename_pop_par, directory, task_name);
% 
% % Analyse the results of the test simulations by finding the accuracy and
% % plotting heatmaps, timecourse predictions and (maybe) RDFs
% results_path = 'workflow_modules/ModelValidation/Results';
% directory = ['workflow_modules/ModelValidation/Results/Figures'];
% analyse_ModelPredictions(obs_range_tol,frac_within_tol,thresh_acceptSAM,thresh_acceptVarSAM,filename_test_data_summary,filename_pop_par,directory, results_path)
% 
% %% 5. Simulating and plotting scenarios
% task_name = 'scenario';
% results_path = 'workflow_modules/ModelPredictions/Results';
% directory = ['workflow_modules/ModelPredictions/Results'];
% simulate_model_scenarios(filename_test_data_summary, filename_patinfo, filename_pop_par, directory,task_name)
% 
% path_scenario = 'workflow_modules/ModelPredictions/Results/scenario';
% directory = 'workflow_modules/ModelPredictions/Results/Figures';
% 
% plot_model_scenarios(pats_test,path_scenario,directory);
