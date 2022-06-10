%% Calculating the LOCF (last observed carried forward) accuracy using the 
% pre and post RDFs as an input to the SAM and VarSAM

%% 1. Load patient data (RDFs)

current_pats_test = load('\final_tiles\PatientData\pat_summary_test.mat')
current_pats_train = load('\final_tiles\PatientData\pat_summary_training.mat')

Lnorm_obs_pre_all = cat(2,current_pats_test.Lnorm_obs_pre,current_pats_train.Lnorm_obs_pre);
Lnorm_obs_post_all = cat(2,current_pats_test.Lnorm_obs_post,current_pats_train.Lnorm_obs_post);

%% 
SAM_prepost = nan(size(Lnorm_obs_post_all));
VarSAM_prepost = nan(size(Lnorm_obs_post_all));

for pat = 1:length(Lnorm_obs_pre_all)
        
    % calculate SAM
   [SAM, varSAM] = computeSAM_v2(Lnorm_obs_pre_all{pat}, Lnorm_obs_post_all{pat},[],[]);
    % store in vector
    SAM_prepost(pat) = SAM;
    VarSAM_prepost(pat) = varSAM;
    
end

% Get the number of accepted patients
accepted = sum(SAM_prepost>0.7 & (VarSAM_prepost>0.3 & VarSAM_prepost<1/0.3))

percent_accepted = accepted/length(Lnorm_obs_pre_all)
