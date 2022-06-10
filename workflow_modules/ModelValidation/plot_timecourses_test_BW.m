%% Oliver Grimm and Lucy Hutchinson 2021
%
function [h]=plot_timecourses_test_BW(currpat, data_summary,p, q,n_CD8_observed_all, n_CD8_observed_PRE_all, set_ind, pat, directory,h)

CD8_timecourse_mean_curr = nanmean(currpat.CD8_timecourse_all,3);
        CD8_timecourse_mean_curr=squeeze(CD8_timecourse_mean_curr);
        % Stack the tiles and parsets on top of each other to reach a 2d matrix
        temp = CD8_timecourse_mean_curr(:,:,1);
        
        for ll = 2:size(CD8_timecourse_mean_curr,3)
            temp = vertcat(temp,CD8_timecourse_mean_curr(:,:,ll));
        end
        
        CD8_timecourse_all_curr  = temp;
        
        figure
        
        plot(nanmedian(CD8_timecourse_all_curr),'k','linewidth',2)
        
        % Use quantiles and fill
        hold on
        quantiles_curr25 = quantile(CD8_timecourse_all_curr,[0.25,0.75]);
        quantiles_curr05 = quantile(CD8_timecourse_all_curr,[0.05,0.95]);
        
        % Put into a long vector for fill
        if sum(sum(isnan(quantiles_curr25)))
            quantiles25_fill=[quantiles_curr25(1,~isnan(quantiles_curr25(1,:))),fliplr(quantiles_curr25(2,~isnan(quantiles_curr25(2,:))))];
        else
            quantiles25_fill = [quantiles_curr25(1,:),fliplr(quantiles_curr25(2,:))];
        end
        
        if sum(sum(isnan(quantiles_curr05)))
            quantiles05_fill=[quantiles_curr05(1,~isnan(quantiles_curr05(1,:))),fliplr(quantiles_curr05(2,~isnan(quantiles_curr05(2,:))))];
        else
            quantiles05_fill = [quantiles_curr05(1,:),fliplr(quantiles_curr05(2,:))];
        end
        
        fill([1:length(quantiles25_fill)/2,length(quantiles25_fill)/2:-1:1]',quantiles25_fill','k','facealpha',0.2,'linestyle','none')
        fill([1:length(quantiles05_fill)/2,length(quantiles05_fill)/2:-1:1]',quantiles05_fill','k','facealpha',0.2,'linestyle','none')
        
        % Add stars for observed data
        plot(ones(size(n_CD8_observed_all{q}))*data_summary.endTimes(q),n_CD8_observed_all{q},'k*')
        % Add CD8S at the start as well
        plot(zeros(size(n_CD8_observed_PRE_all{q})),n_CD8_observed_PRE_all{q},'k*');%b*')
        title(['Pat ',num2str(pat)])
        xlabel('time (days)')
        ylabel('CD8 number per tile')
        xlim([0,data_summary.endTimes(q)+10])
        set(gca,'fontsize',12)
        
        if set_ind == 1
            tag = '_test_pop_pars_';
            title_tag =' TEST ';
        else
            tag = ['_test_control_pars_',num2str(set_ind-1)];
            title_tag=[' Control ',num2str(set_ind-1)];
        end
        
        title(['Patient ',num2str(pat),', CD8 timecourse predictions'],['Parameter set ',title_tag], 'FontSize',14)
        mkdir(directory)
        print([directory,'/TimecourseCD8_simple_median_iqr',tag,'_pat',num2str(pat),'largefontBW'],'-dpng')
        
   %%
%    % Create layered figure
%    if set_ind == 2
%             figure
%             % Plot the timecourse with 25-75% quantiles
%             fill([1:length(quantiles25_fill)/2,length(quantiles25_fill)/2:-1:1]',quantiles25_fill','k','facealpha',0.2,'linestyle','none')
%             hold on
%             plot(nanmedian(CD8_timecourse_all_curr),'k','linewidth',1)
%             h{p}= gca;
%             
%         elseif set_ind ==1
%             
%             fill(h{p},[1:length(quantiles25_fill)/2,length(quantiles25_fill)/2:-1:1]',quantiles25_fill','r','facealpha',0.2,'linestyle','none')
%             hold on
%             plot(h{p},nanmedian(CD8_timecourse_all_curr),'r','linewidth',2)
%             
%             % Add blue stars
%             plot(ones(size(n_CD8_observed_all{q}))*data_summary.endTimes(q),n_CD8_observed_all{q},'k*')
%             plot(zeros(size(n_CD8_observed_PRE_all{q})),n_CD8_observed_PRE_all{q},'k*')
%         else
%             fill(h{p},[1:length(quantiles25_fill)/2,length(quantiles25_fill)/2:-1:1]',quantiles25_fill','k','facealpha',0.2,'linestyle','none')
%             hold on
%             plot(nanmedian(CD8_timecourse_all_curr),'k','linewidth',1)
%             
%    end

