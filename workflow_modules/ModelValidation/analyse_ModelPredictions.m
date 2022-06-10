%% Oliver Grimm and Lucy Hutchinson 2021
%
function analyse_ModelPredictions(obs_range_tol,frac_within_tol,thresh_acceptSAM,...
    thresh_acceptVarSAM,filename_test_datasummary,filename_pop_pars,directory, results_path)

%% Analyse Model Predictions
% Read in the results of the model predictions, calculate the accuracy of
% the predictions and plot the results

%% Load the patient data including post-treat RDFs
%% 1. Get patient details, RDFs and endtimes
load(filename_test_datasummary)
pats_chosen = data_summary.PatID;
n_pats= length(pats_chosen);

% Load the parameter set used for simulations
load(filename_pop_pars);

% Set up variables
n_pars = size(pop_param_vals,1); 
nCD8_endtime_allpats = cell(1,length(pats_chosen));

% If thresholds are not provided, use default
if isempty(obs_range_tol)
    obs_range_tol = 0.2;
end
if isempty(frac_within_tol)
    frac_within_tol = 0.8;
end

total_sets = size(control_param_vals,3)+1;

Lnorm_sim_all = cell(total_sets); % Set up a cell of cells

SAM_matrix = nan(n_pars,length(pats_chosen),total_sets); % Matrix to collect the acceptance of par sets for SAM
VarSAM_matrix = nan(n_pars,length(pats_chosen),total_sets); % Matrix to collect the acceptance of par sets for VarSAM based on range of first 15 distance units
SAM_matrix_NANIndicator = nan(n_pars,length(pats_chosen),total_sets);
mean_completedruns_flags = nan(n_pars,length(pats_chosen),total_sets); % collect run flags (will be 0 if not completed, 1 if completed to endtime)

%% 2. Calculate the RDF and SAM for each patient for each par set
tasknames = {'test', 'control1','control2','control3','control4','control5'};

for set_ind = [2,3,4,5,6,1]
    
    Lnorm_sim_all{set_ind}=cell(length(pats_chosen)); % Save all RDFs from simulations in a cell
    
    for p =1:length(pats_chosen)
        pat = pats_chosen(p);
        q = find(data_summary.PatID == pat); % Find the index of the observed data that corresponds to 'pat'
        
        disp(['SAM loop: patient  ',num2str(pat),' (check pat = ',num2str(data_summary.PatID(q)),')'])
        currpat = load([results_path,'/',tasknames{set_ind},'_pat',num2str(pat)]);
        
        % Calculate the RDF for simulated tiles (all stoch runs)
        n_tiles_curr = length(currpat.tile_inds);
        n_stoch = size(currpat.CD8_timecourse_all,3);
        
        % Get the RDFs of simulated results
        Lnorm_sim_currpat = nan(50, n_tiles_curr,n_stoch,n_pars);
        parfor jj = 1:n_pars % for each parameter set
            kk=0;
            for i_tile = 1:n_tiles_curr
                for i_stoch = 1:n_stoch
                    kk = kk+1;
                    grid_curr = currpat.CD8gridstack_endtime_all{1, jj}{i_tile, i_stoch};
                    Lnorm_sim_currpat(:,i_tile,i_stoch,jj) = computeRDF(grid_curr,50);
                end
            end
        end
        
        % Reshape Lnorm_sim
        Lnorm_sim_currpat_reshape=Lnorm_sim_currpat(:,:,1,:);
        for i_stoch = 2:n_stoch
            Lnorm_sim_currpat_reshape=[Lnorm_sim_currpat_reshape,Lnorm_sim_currpat(:,:,i_stoch,:)];
        end
        Lnorm_sim_all{set_ind}{p} = squeeze(Lnorm_sim_currpat_reshape);
        
        % Calculate the SAM and store in matrix
        Lnorm_obs_post_curr = Lnorm_obs_post{q};
        parfor jj = 1:n_pars % for each parameter set
            mean_completedruns_flags(jj,p,set_ind) = mean(mean(currpat.flags_finished_all{jj}));% save the completed run flags
            Lnorm_sim_currpar = Lnorm_sim_all{set_ind}{p}(:,:,jj);
            [SAM_out,VarSAM_out] = computeSAM(Lnorm_sim_currpar, Lnorm_obs_post_curr, obs_range_tol, frac_within_tol);
            SAM_matrix(jj,p, set_ind) = SAM_out;
            VarSAM_matrix(jj,p,set_ind) = VarSAM_out;
        end
        
    end % pats
    
    % Replace the entries in the acceptance matrices with NaNs if there were
    % not enough non-NaN results to compare the RDFs.
    maxNaNs_frac = 0.5; % Maximum fraction of NaNs (resulting from too few CD8s) that means SAM can't be calculated
    
    [SAM_matrix_NANIndicator_curr]=indicate_nans_SAM(pats_chosen,data_summary,Lnorm_obs_post,Lnorm_sim_all{set_ind},SAM_matrix(:,:,set_ind),maxNaNs_frac);
    SAM_matrix_NANIndicator(:,:,set_ind)=SAM_matrix_NANIndicator_curr;
    
end   % sets

%% Plot barplots and get a value for the accuracy
% To get a numerical value for the success of the test, find the mean SAM
% across all tested pars

% Make VarSAM<1

VarSAM_scale = VarSAM_matrix;
VarSAM_scale(VarSAM_scale>1) = 1./(VarSAM_scale(VarSAM_scale>1)); % for values >1 (sim range greater than obs range), flip the fraction so it is all on the same scale


meanSAM1 = nanmean(nanmean(SAM_matrix_NANIndicator,2),1);
meanVarSAM = nanmean(nanmean(VarSAM_scale,2),1);

SAM_votes_matrix =(SAM_matrix_NANIndicator>=thresh_acceptSAM).*(VarSAM_scale>=thresh_acceptVarSAM);

test_results_condensed = reshape(mean(SAM_votes_matrix),n_pats,[]);

save([directory,'/processed_testresults'],'test_results_condensed','SAM_votes_matrix','VarSAM_scale','meanVarSAM','meanSAM1',...
    'pats_chosen','SAM_matrix_NANIndicator','VarSAM_matrix','Lnorm_sim_all');

% Make a bar plot
plot_barplot(test_results_condensed, pats_chosen, total_sets,[directory,'/barplot_test'])

%% Plot Heatmaps
plot_SAM_heatmaps_test(pats_chosen, SAM_matrix_NANIndicator,VarSAM_matrix, pop_param_vals, par_names, total_sets,[directory,'/SAM_heatmap_test'])


%% Now produce plots of timecourses
h={};
for set_ind =1% [2,3,4,5,6,1]%1:total_sets
    
    
    for p =1:length(pats_chosen)
        pat = pats_chosen(p);
        q = find(data_summary.PatID == pat); % Find the index of the observed data that corresponds to 'pat'
        currpat = load([results_path,'/',tasknames{set_ind},'_pat',num2str(pat)]);
        h=plot_timecourses_test_BW(currpat, data_summary,p, q,n_CD8_observed_all, n_CD8_observed_PRE_all, set_ind, pat, directory,h);
    end
end


