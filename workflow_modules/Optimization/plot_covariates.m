%% Oliver Grimm and Lucy Hutchinson 2021
%
function []= plot_covariates(pats_chosen, data_summary_in, mean_completedruns_flags, SAM_matrix_NANIndicator,VarSAM_scale,thresh_acceptSAM,thresh_acceptVarSAM, pars_all, filename)

%% A function to plot covariates against 
% 1. fraction of completed runs for each pat
% 2. Fraction of accepted par sets according to SAM

n_pars = length(pars_all);
data_summary = data_summary_in(ismember(data_summary_in.PatID,pats_chosen),:);

%% 1. Use the mean flags to get a fraction of completed runs (across all tested par sets) for each patient
fractionCompleteruns = [sum(mean_completedruns_flags==1)/n_pars]';



%% 2. Use SAM to get fraction of accepted runs according to SAM (at some threshold thresh) 
if isempty(thresh_acceptSAM)
    thresh_acceptSAM = 0.7;
end
if isempty(thresh_acceptVarSAM)
    thresh_acceptVarSAM = 0.2;
end

acceptSAM = SAM_matrix_NANIndicator>thresh_acceptSAM & VarSAM_scale >thresh_acceptVarSAM;
fractionAcceptSAM = sum(logical(acceptSAM))'/n_pars;


%% Plot all figures in the same way
for plot_ID = 1:2
    if plot_ID ==1
        xdata = fractionCompleteruns;
        xdata_name = 'fraction completed runs';
        print_name = 'fracCompleteRuns';
    elseif plot_ID ==2
        xdata = fractionAcceptSAM;
        xdata_name = 'fraction accepted pars SAM'; 
        print_name = 'fracAcceptedParsSAM';
    end

figure('units','normalized','position',[0.02 0.02 0.95 0.95])
subplot_ind=0;
for ii = [2:7,9:10,12:14]
    subplot_ind = subplot_ind+1;
    subplot(3,4,subplot_ind)
    plot(xdata,data_summary.(ii),'o')
    hold on 
    lsline
    xlabel(xdata_name)
    ylabel(data_summary.Properties.VariableNames{ii})  
end
print([filename,'covariates_',print_name],'-dpng')

end
