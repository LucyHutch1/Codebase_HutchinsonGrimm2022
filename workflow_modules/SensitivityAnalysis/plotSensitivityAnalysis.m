%% Oliver Grimm and Lucy Hutchinson 2021
%
%% Plotting the results of the sensitivity analysis
function plotSensitivityAnalysis(filename_in, directory,pats_chosen, tiles_chosen)
% 1. Read in the results
load(filename_in);

% Depending on how many parameters were scanned, make several figures. Each
% figure contains results from 4 parameters
num_parsperfig = 4;
numfigs = ceil(size(SA_allData,1)/num_parsperfig);

filenamef1=cell(numfigs);
filenamef2=cell(numfigs);
filenamef3=cell(numfigs);
filenamef4=cell(numfigs);

mkdir(directory)

for ind_filename = 1:numfigs
    
    filenamef1{ind_filename} = ([directory,'/CellNumber_time_',num2str(ind_filename)]);
    filenamef2{ind_filename} = ([directory,'/RDF_',num2str(ind_filename)]);
    filenamef3{ind_filename} = ([directory,'/CellNumber_vs_ParVal_t20_',num2str(ind_filename)]);
    filenamef4{ind_filename} = ([directory,'/CellNumber_vs_ParVal_t50_',num2str(ind_filename)]);
    
end

% Set up variables
CD8curr_pat = zeros(size(SA_allData,2), size(SA_allData{1,1,1}{4},2));
TUMcurr_pat = zeros(size(SA_allData,2), size(SA_allData{1,1,1}{4},2));

par_vals_curr= zeros(size(SA_allData,2),1);
CD8curr_pat_last=zeros(size(SA_allData,2),1);
TUMcurr_pat_last=zeros(size(SA_allData,2),1);
Lcurr_pat = zeros(size(SA_allData,2), size(SA_allData{1,1,1}{7},2)); % L corresponds to RDF

% One column for each starting tile
% One row for each varied paramter
% Colours to represent values of the parameter

num_tiles = size(SA_allData,3);
num_pars = size(SA_allData,1);

plot_ind = 0;

%% f1. CD8 number over time

f1{1}=figure('units','centimeters','position',[0 0 35 25]);
ind_curr_fig = 1;
for ii = 1:num_pars % loop over parameters
    for jj = 1:num_tiles % loop over tiles
        % Put simulations for different par values in a matrix
        for kk = 1:size(SA_allData,2)
            
            %% Use median instead of mean to take into account bistability
            CD8curr_pat(kk,:) = nanmedian(SA_allData{ii,kk,jj}{8},1);
            TUMcurr_pat(kk,:) = nanmedian(SA_allData{ii,kk,jj}{9},1);
            
        end
        plot_ind = plot_ind+1;
        if plot_ind > num_parsperfig*num_tiles
            plot_ind = 1;
            % Add title including which page we are up to
            annotation('textbox',[0.38,0.97,0.1,0.02],'String',['Sensitivity analysis: median cell number over time (',num2str(ind_curr_fig),'/',num2str(numfigs),')'],'FitBoxToText','on');
            ind_curr_fig = ind_curr_fig +1;
            f1{ind_curr_fig}=figure('units','centimeters','position',[0 0 35 25]);
        end
        subplot(num_parsperfig,num_tiles,plot_ind)
        % Red shades for tumour cells
        red_cols = autumn(size(TUMcurr_pat,1));
        % Blue shades for CD8s
        blue_cols = winter(size(CD8curr_pat,1));
        for ll = 1:size(CD8curr_pat,1)
            plot(1:size(CD8curr_pat,2), CD8curr_pat(ll,:), 'color',blue_cols(ll,:),'linewidth',1.5)
            hold on
            plot(1:size(TUMcurr_pat,2), TUMcurr_pat(ll,:), 'color',red_cols(ll,:),'linewidth',1.5)
        end
        title([SA_allData{ii,1,jj}{1},' Tile ',num2str(jj)])
        xlabel('time')
        ylabel('Cell number')
        xlim([0 size(TUMcurr_pat,2)])
    end
end
% Add annotation to the last figure
annotation('textbox',[0.38,0.97,0.1,0.02],'String',['Sensitivity analysis: median cell number over time (',num2str(ind_curr_fig),'/',num2str(numfigs),')'],'FitBoxToText','on');


%% f3 and f4. Cell numbers versus parameter value time 20 and time 50 (or last observed if <20 or <50)
tpts = [20,50];
for tpt_ind = 1:2
    tpt_curr = tpts(tpt_ind);
    if tpt_ind ==1
        f3{1} = figure('units','centimeters','position',[0 0 35 25]);
    else
        f4{1} = figure('units','centimeters','position',[0 0 35 25]);
        
    end
    plot_ind=0;
    ind_curr_fig = 1;
    for ii = 1:num_pars % loop over parameters
        for jj = 1:num_tiles % loop over tiles
            
            % Put simulations for different par values in a matrix
            for kk = 1:size(SA_allData,2)
                last_tpt = max(find(~isnan(SA_allData{ii,kk,jj}{4})));
                tpt_to_plot = min(tpt_curr,last_tpt);
                CD8curr_pat(kk,:) = SA_allData{ii,kk,jj}{4};
                TUMcurr_pat(kk,:) = SA_allData{ii,kk,jj}{5};
                
                par_vals_curr(kk) = SA_allData{ii,kk,jj}{2};
                CD8curr_pat_last(kk)=CD8curr_pat(kk,tpt_to_plot);
                TUMcurr_pat_last(kk)=TUMcurr_pat(kk,tpt_to_plot);
            end
            
            plot_ind = plot_ind+1;
            
            if plot_ind > num_parsperfig*num_tiles
                plot_ind = 1;
                ind_curr_fig = ind_curr_fig +1;
                if tpt_ind ==1
                    annotation('textbox',[0.37,0.97,0.1,0.02],'String',['Sensitivity analysis: median cell number at 20 days (',num2str(ind_curr_fig-1),'/',num2str(numfigs),')'],'FitBoxToText','on');
                    f3{ind_curr_fig}=figure('units','centimeters','position',[0 0 35 25]);
                else
                    annotation('textbox',[0.37,0.97,0.1,0.02],'String',['Sensitivity analysis: median cell number at 50 days (',num2str(ind_curr_fig-1),'/',num2str(numfigs),')'],'FitBoxToText','on');
                    f4{ind_curr_fig}=figure('units','centimeters','position',[0 0 35 25]);
                end
            end
            subplot(num_parsperfig,num_tiles,plot_ind)
            
            red_cols = autumn(size(TUMcurr_pat,1));
            blue_cols = winter(size(CD8curr_pat,1));
            
            % Plot the parameters that were sampled on a log scale 
            if ii == 7 | ii == 9 | ii==10
            for ll = 1:size(CD8curr_pat,1)
                semilogx(par_vals_curr, CD8curr_pat_last,'o-', 'color',blue_cols(1,:),'linewidth',1.5)
                hold on
                semilogx(par_vals_curr, TUMcurr_pat_last,'o-', 'color',red_cols(1,:),'linewidth',1.5)
            end    
                
                
            else   
            for ll = 1:size(CD8curr_pat,1)
                plot(par_vals_curr, CD8curr_pat_last,'o-', 'color',blue_cols(1,:),'linewidth',1.5)
                hold on
                plot(par_vals_curr, TUMcurr_pat_last,'o-', 'color',red_cols(1,:),'linewidth',1.5)
            end
            end
            title([SA_allData{ii,1,jj}{1},' Tile ',num2str(jj)]);
            xlabel('parameter value')
            ylabel('Cell number')
        end
    end
    
    if tpt_ind ==1
        annotation('textbox',[0.37,0.97,0.1,0.02],'String',['Sensitivity analysis: median cell number at 20 days or last obs (',num2str(ind_curr_fig),'/',num2str(numfigs),')'],'FitBoxToText','on');
    else
        annotation('textbox',[0.37,0.97,0.1,0.02],'String',['Sensitivity analysis: median cell number at 50 days or last obs (',num2str(ind_curr_fig),'/',num2str(numfigs),')'],'FitBoxToText','on');
    end
end


%% f2. Plotting RDF at the last observed timepoint
plot_ind = 0;
ind_curr_fig = 1;
f2{1} = figure('units','centimeters','position',[0 0 35 25]);
cool_cols = cool(size(Lcurr_pat,1));
for ii = 1:num_pars % loop over parameters
    for jj = 1:num_tiles % loop over tiles
        % Put simulations for different par values in a matrix
        for kk = 1:size(SA_allData,2)
            Lcurr_pat(kk,:) = SA_allData{ii,kk,jj}{7};
        end
        plot_ind = plot_ind+1;
        if plot_ind > num_parsperfig*num_tiles
            plot_ind = 1;
            annotation('textbox',[0.38,0.97,0.1,0.02],'String',['Sensitivity analysis: RDF at 100 days (',num2str(ind_curr_fig),'/',num2str(numfigs),')'],'FitBoxToText','on');
            ind_curr_fig = ind_curr_fig +1;
            f2{ind_curr_fig}=figure('units','centimeters','position',[0 0 35 25]);
        end
        subplot(num_parsperfig,num_tiles,plot_ind)
        
        for ll = 1:size(Lcurr_pat,1)
            plot(1:50, Lcurr_pat(ll,:), 'color',cool_cols(ll,:),'linewidth',1.5)
            hold on
        end
        title([SA_allData{ii,1,jj}{1},' Tile ',num2str(jj)])
        xlabel('distance')
        ylabel('Rel cell density')
        xlim([0 50])
        
    end
end
annotation('textbox',[0.38,0.97,0.1,0.02],'String',['Sensitivity analysis: RDF at 100 days (',num2str(ind_curr_fig),'/',num2str(numfigs),')'],'FitBoxToText','on');

%% Save figures
for ind_fig = 1:numfigs
    print(f1{ind_fig},filenamef1{ind_fig},'-dpng')
    print(f2{ind_fig},filenamef2{ind_fig},'-dpng') % f2 is not used in the manuscript
    print(f3{ind_fig},filenamef3{ind_fig},'-dpng') 
    print(f4{ind_fig},filenamef4{ind_fig},'-dpng')
end



%%%%%%%%%%%%%%
%% Bar Plot %%
%%%%%%%%%%%%%%

% Tabulate baseline numbers of CD8 and Tum cells (need to load actual
% structs as we don't have number of tumour cells in a table)
CD8_baseline = nan(1,6); TUM_baseline = nan(1,6);
for p = 1:6
    pat = pats_chosen(p);
    % read in pre- structs
        filename_preStruct = ['final_tiles/PatientData/p',num2str(pat),'/','pre','_',...
            '10','percent/structs_p',num2str(pat),'pre.mat']; 
    load(filename_preStruct);
    curr_tile = allSystems{tiles_chosen(p)};
    CD8_baseline(p) = sum(sum(curr_tile.grid.Li));
    TUM_baseline(p) = sum(sum(curr_tile.grid.Lt));
end

% Calculate delta(CD8)/delta(Par) and delta(TUM)/delta(Par) as the
% sensitivity measure

CD8_sensitivity = nan(13,6);
TUM_sensitivity = nan(13,6);

% Look at the 20 day timepoint
querytime = 20;

% Loop through SA results to prepare the sensitivity measures of each
% parameter
for j = 1:num_pars
    for p = 1:6
        % We wish to find the parameter range scanned so subtract the
        % smallest value for each parameter (1) from the largest (10)
        curr_res1 = SA_allData{j,1,p};
        curr_res10 = SA_allData{j,10,p};
        
        % Get last non NAN timepoint
        last_tpt1_CD8(j,p) = max(find(~isnan(curr_res1{4})));
        last_tpt10_CD8(j,p) = max(find(~isnan(curr_res10{4})));
        last_tpt1_TUM(j,p) = max(find(~isnan(curr_res1{5})));
        last_tpt10_TUM(j,p) = max(find(~isnan(curr_res10{5})));
        
        delta_CD8(j,p) = curr_res1{4}(min(last_tpt1_CD8(j,p),querytime)) - curr_res10{4}(min(last_tpt10_CD8(j,p),querytime)); % 4th entry is the CD8 number tile mean at all time points
        delta_TUM(j,p) = curr_res1{5}(min(last_tpt1_TUM(j,p),querytime)) - curr_res10{5}(min(last_tpt10_TUM(j,p),querytime)); % 5th entry is tumor tile mean
        cfb_CD8(j,p) = delta_CD8(j,p)/CD8_baseline(p); % change from baseline
        cfb_TUM(j,p) = delta_TUM(j,p)/TUM_baseline(p); % change from baseline
        
        delta_par = curr_res10{2}- curr_res1{2}; 
        midpt_par = median([curr_res10{2}, curr_res1{2}]); % midpoint to normalise by
        
        CD8_sensitivity(j,p)= cfb_CD8(j,p)/(delta_par/midpt_par);
        TUM_sensitivity(j,p) = cfb_TUM(j,p)/(delta_par/midpt_par);
    end
pars{j}=SA_allData{j,1,1}{1};
end

%

figure('units','normalized','position',[0.1,0.1,0.8,0.5])
subplot(1,2,1)
b1 = bar(log(abs(CD8_sensitivity)+1))
xticklabels(pars)
xtickangle(45)
ylabel('CD8 sensitivity')
set(gca,'fontsize',14)
getylims = get(gca,'ylim')
title(['Sensitivity of CD8 cell number at ',num2str(querytime),' days'])

subplot(1,2,2)
bar(log(abs(TUM_sensitivity)+1))
xticklabels(pars)
xtickangle(45)
ylabel('Tumour sensitivity')
set(gca,'fontsize',14)
ylim(getylims)
title(['Sensitivity of tumour cell number at ',num2str(querytime),' days'])
legend({['CD8 baseline = ',num2str(CD8_baseline(1))],['CD8 baseline = ',num2str(CD8_baseline(2))],...
    ['CD8 baseline = ',num2str(CD8_baseline(3))],['CD8 baseline = ',num2str(CD8_baseline(4))],...
    ['CD8 baseline = ',num2str(CD8_baseline(5))],['CD8 baseline = ',num2str(CD8_baseline(6))]}, 'location','northeast')


print([directory,'/Barplot_',num2str(querytime),'days'],'-dpng')