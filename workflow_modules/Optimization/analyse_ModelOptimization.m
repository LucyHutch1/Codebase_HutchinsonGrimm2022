%% Oliver Grimm and Lucy Hutchinson 2021
%
function []=analyse_ModelOptimization(obs_range_tol, frac_within_tol,datasummary_path,pats_chosen,parsamples_path, results_path,task_name, savepath)
% Processing model optimization results
% Read in results from the optimization and calculate the SAM for each
% parameter set for each patient

%% 1. Get patient details, RDFs and endtimes
load(datasummary_path);
pats = pats_chosen;

% Load the parameter set used for simulations
load(parsamples_path);

% Set up variables
n_pars = size(pars_all,1); % Don't hard code 
nCD8_endtime_allpats = cell(1,length(pats_chosen));

if isempty(obs_range_tol)
    obs_range_tol = 0.2;
end
if isempty(frac_within_tol)
    frac_within_tol = 0.8;
end

Lnorm_sim_all = cell(1,length(pats_chosen)); % Save all RDFs from simulations in a cell
SAM_matrix = nan(n_pars,length(pats_chosen)); % Matrix to collect the acceptance of par sets for SAM
VarSAM_matrix = nan(n_pars,length(pats_chosen)); % Matrix to collect the acceptance of par sets for VarSAM based on range of first 15 distance units
mean_completedruns_flags = nan(n_pars,length(pats_chosen)); % collect run flags (will be 0 if not completed, 1 if completed to endtime)

%% 2. Calculate the RDF and SAM for each patient for each par set
for p =1:length(pats_chosen)
    pat = pats_chosen(p);
    q = find(data_summary.PatID == pat); % Find the index of the observed data that corresponds to 'pat'
    
    disp(['SAM loop: patient  ',num2str(pat),' (check pat = ',num2str(data_summary.PatID(q)),')'])
    currpat = load([results_path,'/',task_name,'_pat',num2str(pat)]);

    % Calculate the RDF for simulated tiles (all stoch runs)
    n_tiles_curr = length(currpat.tile_inds);
    n_stoch = 5;
    
    % Get the RDFs of simulated results %PARALLELIZE HERE
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
        Lnorm_sim_all{p} = squeeze(Lnorm_sim_currpat_reshape);
        Lnorm_sim_currpat_reshape = squeeze(Lnorm_sim_currpat_reshape); 
        
    % Calculate the SAM and store in matrix
    Lnorm_obs_post_curr = Lnorm_obs_post{q};
    parfor jj = 1:n_pars % for each parameter set
        mean_completedruns_flags(jj,p) = mean(mean(currpat.flags_finished_all{jj}));% save the completed run flags
        Lnorm_sim_currpar = Lnorm_sim_currpat_reshape(:,:,jj);
        [SAM_out,VarSAM_out] = computeSAM(Lnorm_sim_currpar, Lnorm_obs_post_curr, obs_range_tol, frac_within_tol);
        SAM_matrix(jj,p) = SAM_out;
        VarSAM_matrix(jj,p) = VarSAM_out;
    end
end % patient loop

% Replace the entries in the acceptance matrices with NaNs if there were
% not enough non-NaN results to compare the RDFs.
maxNaNs_frac = 0.5; % Maximum fraction of NaNs (resulting from too few CD8s) that means SAM can't be calculated

[SAM_matrix_NANIndicator]=indicate_nans_SAM(pats_chosen,data_summary,Lnorm_obs_post,Lnorm_sim_all,SAM_matrix,maxNaNs_frac);

% Scale VarSAM so it falls between zero and one
VarSAM_scale = VarSAM_matrix;
VarSAM_scale(VarSAM_scale>1) = 1./(VarSAM_scale(VarSAM_scale>1)); % for values >1 (sim range greater than obs range), flip the fraction so it is all on the same scale



%% 3. Save KS and SAM tables that can be used as an input to the virtual population function
 save([savepath,'.mat'],'SAM_matrix','VarSAM_matrix','VarSAM_scale', 'Lnorm_sim_all','Lnorm_obs_post','data_summary','pats_chosen','SAM_matrix_NANIndicator','mean_completedruns_flags');
 % Also save a light version that doesn't include the RDFs (easier to
 % store and load)
 save([savepath,'_lite.mat'],'SAM_matrix','VarSAM_matrix','VarSAM_scale','data_summary','pats_chosen','SAM_matrix_NANIndicator');
