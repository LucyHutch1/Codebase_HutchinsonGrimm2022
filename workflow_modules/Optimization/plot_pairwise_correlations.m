%% Oliver Grimm and Lucy Hutchinson 2021
%
function [] = plot_pairwise_correlations(pats_chosen,pars_all, par_names,SAM_matrix,VarSAM_scale, thresh_acceptSAM, thresh_acceptVarSAM, thresh_votesforaccept, save_path)
%% function to plot the parameters in a pairwise fashion and colour by the votes from SAM

% If no threshold is given
if isempty(thresh_acceptSAM)
    thresh_acceptSAM = 0.7;
end
if isempty(thresh_acceptVarSAM)
    thresh_acceptVarSAM = 0.2;
end


% Colour by SAM
% Make the SAM matrix into a logical by accepting any values over a
% threshold for SAM and VarSAM
acceptSAM = SAM_matrix>thresh_acceptSAM & VarSAM_scale > thresh_acceptVarSAM;
sum_votes = sum(logical(acceptSAM),2);
pop_param_ind = find(sum_votes>thresh_votesforaccept*length(pats_chosen));
    
title_str = ['Coloured by SAM votes (thresh = [',num2str(thresh_acceptSAM),',',num2str(thresh_acceptVarSAM),'])'];

    
figure('units','normalized','position',[0 0 1 1])
par_pairs = [1,2;1,3;1,4;2,3;2,4;3,4];
for p = 1:size(pats_chosen,2)
    for kk = 1:6
        subplot(2,3,kk)
        if par_pairs(kk,2) == 4
        scatter(log10(pars_all(:,par_pairs(kk,1))), pars_all(:,par_pairs(kk,2)),[],sum_votes)
        hold on
        scatter(log10(pars_all(pop_param_ind,par_pairs(kk,1))), pars_all(pop_param_ind,par_pairs(kk,2)),[],sum_votes(pop_param_ind),'filled','MarkerEdgeColor','r')
        xlabel(['log10(',par_names{par_pairs(kk,1)},')'])
        ylabel(par_names{par_pairs(kk,2)})
        ax=gca;
        ax.FontSize=14;
        else
        scatter(log10(pars_all(:,par_pairs(kk,1))), log10(pars_all(:,par_pairs(kk,2))),[],sum_votes)
        hold on
        scatter(log10(pars_all(pop_param_ind,par_pairs(kk,1))), log10(pars_all(pop_param_ind,par_pairs(kk,2))),[],sum_votes(pop_param_ind),'filled','MarkerEdgeColor','r')
        xlabel(['log10(',par_names{par_pairs(kk,1)},')'])
        ylabel(['log10(',par_names{par_pairs(kk,2)},')'])
        ax=gca;
        ax.FontSize=14;  
        end
        
    end

end

c=colorbar('manual','position',[0.93,0.3,0.02,0.4]);
c.Label.String='Votes';
c.Label.FontSize=15;
print([save_path,'/scatterplots_SAMvotes'],'-dpng')

