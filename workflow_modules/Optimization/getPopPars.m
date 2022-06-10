%% Oliver Grimm and Lucy Hutchinson 2021
%
function getPopPars(thresh_acceptSAM, thresh_acceptVarSAM,thresh_votesforaccept,opt_results_processed, filename_parsets,filename_datasummary, results_path, task_name)

% Using the optimization results, get the best performing parameter sets
load([opt_results_processed,'_lite']);
load(filename_parsets);
load(filename_datasummary)
n_pars = length(pars_all);

% If there are no thresholds provided, use the defaults
if isempty(thresh_acceptSAM)
thresh_acceptSAM = 0.7;
end
if isempty(thresh_acceptVarSAM)
thresh_acceptVarSAM = 0.3;
end
if isempty(thresh_votesforaccept)
thresh_votesforaccept = 0.6;
end


% Get all the 'definitely accepted' pop parameters by looking at the votes
% for each parameter set 
acceptSAM = SAM_matrix_NANIndicator>thresh_acceptSAM & VarSAM_scale > thresh_acceptVarSAM;

%% Rescue 
% For patients who had too few CD8s in on-treatment samples, use
% the number of CD8 cells to find the best performing parameter sets. 
% Find patients for which this is the case

NaNPats = pats_chosen(sum(isnan(SAM_matrix_NANIndicator))==n_pars);

for ii = 1:length(NaNPats)
    pat = NaNPats(ii);
    q = find(data_summary.PatID==pat);
    p = find(pats_chosen==pat);
    CD8_obs = n_CD8_observed_all{q};
    currpat_optresults=load([results_path,'/',task_name,'_pat',num2str(pat)]);
    for parset = 1:n_pars
    CD8_sim = reshape(currpat_optresults.nCD8_endtime_all{parset},[],1);
    [h,prob]=ttest2(CD8_obs,CD8_sim);
    % Add accept/reject to the acceptSAM matrix in the correct place. H=0
    % means distributions are the same so accept the parameter set.
    acceptSAM(parset,p)=~logical(h);
    end
end

% Also attempt to save individual cases where a given parameter set
% resulted in too few CD8s in the simulated tiles
% Make a copy of the SAM matrix where the NaNPats no longer have NaNs
%%
SAM_matrix_NANIndicator_copy = SAM_matrix_NANIndicator;
SAM_matrix_NANIndicator_copy(:,sum(isnan(SAM_matrix_NANIndicator_copy))==n_pars)=2;
[nanpar_ind,nanpat_ind]=find(isnan(SAM_matrix_NANIndicator_copy));
nanpat_unique=unique(nanpat_ind);

for ii = 1:length(nanpat_unique)
    
    p=nanpat_unique(ii);
    pat = pats_chosen(p);
    q = find(data_summary.PatID==pat);
    CD8_obs = n_CD8_observed_all{q};
    currpat_optresults=load([results_path,'/',task_name,'_pat',num2str(pat)]);
    nanpar_ind_curr = nanpar_ind(nanpat_ind == p);
    for parset = nanpar_ind_curr'
    CD8_sim = reshape(currpat_optresults.nCD8_endtime_all{parset},[],1);
    [h,~]=ttest2(CD8_obs,CD8_sim);
    % Add accept/reject to the acceptSAM matrix in the correct place. H=0
    % means distributions are the same so accept the parameter set.
    acceptSAM(parset,p)=~logical(h);
    end
    
end



%%
sum_votes = sum(logical(acceptSAM),2);
pop_param_ind = find(sum_votes>thresh_votesforaccept*length(pats_chosen));
pop_param_vals = pars_all(pop_param_ind,:);

% Make some randomly selected control parameter sets as well (selected from
% optimisation pars)
control_param_vals = nan(length(pop_param_ind), 4, 5);
for ii = 1:5
    control_param_vals(:,:,ii) = pars_all(randperm(length(pars_all), length(pop_param_ind)), :);
end
 save([results_path,'/pop_params'], 'acceptSAM','sum_votes','pop_param_ind','pop_param_vals', 'par_names','control_param_vals')   