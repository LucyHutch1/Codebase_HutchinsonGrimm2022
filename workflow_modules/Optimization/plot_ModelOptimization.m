%% Oliver Grimm and Lucy Hutchinson 2021
%
function plot_ModelOptimization(thresh_acceptSAM,thresh_acceptVarSAM,thresh_votesforaccept, obs_range_tol, opt_results_processed, par_samples_opt, savepath)% Plot optimization results
load(opt_results_processed);
load(par_samples_opt);

% If the thresholds are not provided, use the defaults
if isempty(thresh_acceptSAM)
thresh_acceptSAM = 0.7; 
end
if isempty(thresh_acceptVarSAM)
    thresh_acceptVarSAM = 0.3; 
end
if isempty(thresh_votesforaccept)
thresh_votesforaccept = 0.6;
end
if isempty(obs_range_tol)
 obs_range_tol = 0.2;
end

% 1. Heatmaps of SAM
filename_heatmaps = [savepath,'heatmaps'];
plot_SAM_heatmaps(pats_chosen, SAM_matrix_NANIndicator, VarSAM_scale, pars_all, par_names, filename_heatmaps);

% 1. Scatter plots for covariates
plot_covariates(pats_chosen, data_summary, mean_completedruns_flags, SAM_matrix_NANIndicator,VarSAM_scale,thresh_acceptSAM,thresh_acceptVarSAM, pars_all, savepath);

% 3. SAM raw data for all par sets (Note: this will generate 20-40 pages
% per patient, so just run for one example)
 p=3;

 plot_RDF_SAM(Lnorm_sim_all{p}, Lnorm_obs_post{p}, obs_range_tol, SAM_matrix_NANIndicator(:,p),VarSAM_matrix(:,p), thresh_acceptSAM,thresh_acceptVarSAM,1, pats_chosen(p),savepath)

% 4. Colourful scatter plot to identify pairwise correlations between
% parameters
plot_pairwise_correlations(pats_chosen,pars_all, par_names,SAM_matrix_NANIndicator,VarSAM_scale, thresh_acceptSAM, thresh_acceptVarSAM, thresh_votesforaccept, savepath)
