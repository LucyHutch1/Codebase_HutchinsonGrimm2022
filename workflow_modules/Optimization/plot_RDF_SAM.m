%% Oliver Grimm and Lucy Hutchinson 2021
%
function [] = plot_RDF_SAM(Lnorm_sim, Lnorm_obs, obs_range_tol,SAM_matrix_NANIndicator,VarSAM_matrix,thresh_acceptSAM, thresh_acceptVarSAM, save_plots, pat, save_path)
 
% Find the max and min and range observed for each value of x (1,50) (same
% as where we calculated the SAM)
max_obs = max(Lnorm_obs,[],2);
min_obs = min(Lnorm_obs,[],2);
range_obs = max_obs-min_obs;

% If no threshold is given
if isempty(thresh_acceptSAM)
    thresh_acceptSAM = 0.7;
end
if isempty(thresh_acceptVarSAM)
    thresh_acceptVarSAM = 0.2;
end


% Get the upper and lower thresholds for acceptance, based on the tolerance
% (default 0.2)
upper_thresh = max_obs+obs_range_tol*range_obs;
lower_thresh = max(0,min_obs-obs_range_tol*range_obs);

% Make a subfolder
mkdir([save_path,'/RDF_SAM_plots_pat',num2str(pat)])

fig_num = 0;
figure('units','normalized','position',[0 0 1 1])
for par_set = 1:size(Lnorm_sim,3)
    subplot_ind = rem(par_set,25);
    
    if subplot_ind == 0
        subplot_ind = 25;
    end
    subplot(5,5,subplot_ind)
    % Plot the bounds with thick blue lines
    plot([1:50],Lnorm_sim(:,:,par_set),'r','linewidth',0.3)
    hold on
    fill([1:50,50:-1:1]', [upper_thresh; flipud(lower_thresh)],'b','facealpha',0.2,'linestyle','none')
    plot([1:50], upper_thresh,'b-','linewidth',1.5)
    plot([1:50], lower_thresh,'b-','linewidth',1.5)
    title(['Param set ',num2str(par_set)])
    
    % Shade the background mint green for accepted parameter sets
    if SAM_matrix_NANIndicator(par_set)>thresh_acceptSAM & VarSAM_matrix >thresh_acceptVarSAM
        currax = gca;
        currax.Color = [225, 255, 243]/255;
    end
    
    
   if subplot_ind == 25 
        if save_plots ==1
            % Save figures only if the switch is on
            fig_num = fig_num+1;

            print([save_path,'/RDF_SAM_plots_pat',num2str(pat),'/RDF_SAM_figs_pat',num2str(pat),'_fig',num2str(fig_num)],'-dpng')

        end
            figure('units','normalized','position',[0 0 1 1])
    end 
end
