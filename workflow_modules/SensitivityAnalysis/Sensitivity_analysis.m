%% Oliver Grimm and Lucy Hutchinson 2021
%
%% Sensitivity analysis
% A function to test the sensitivity of model parameters
function Sensitivity_analysis(n_stoch,pats_chosen, tiles_chosen ,directory,filenameout)

% Set physiological ranges for each of the 13 model parameters and store in
% a cell array with the names of the parameters
pars = {'TUpprol', 'TUpmig', 'TUpdeath', 'TUpmax', 'TUintmax', 'IMpmax', 'IMpprol', 'IMpmig', 'IMpkill', 'IMpdeath', 'IMrwalk', 'IMinfluxProb', 'IMinflRate'};
ranges = [linspace(0.05,0.5,10);... % 1 TUpprol
    linspace(0,0.35,10);... % 2 TUpmig
    linspace(0.01,0.12,10);...% 3 TUpdeath
    6:15;...% 4 TUpmax
    1:10;...% 5 TUintmax
    3:12;...% 6 IMpmax
    logspace(-2,-0.82,10);... % 7 IMpprol
    linspace(0.05,0.7,10);... % 8 IMpmig
    logspace(-3,-0.52,10);... % 9 IMpkill
    logspace(-3,-1.12,10);...% 10 IMpdeath
    linspace(0.5,1,10);... % 11 IMrwalk
    linspace(0.01,0.5,10);... % 12 IMinfluxProb
    1:10]; % 13 IMinflRate

param_array = {};
for ii = 1:length(pars)
   param_array{ii,1} = pars{ii};
   param_array{ii,2} = ranges(ii,:);
end


% Set simulation end time
endTime = 100;

% Set up an empty cell array for the results
SA_allData = {}; 


%% Run simulations for all parameter sets for each patient
for p = 1:length(pats_chosen) % patient loop
    pat = pats_chosen(p);
    
    % read in pre- structs for current patient
    filename_preStruct = ['final_tiles/PatientData/p',num2str(pat),'/','pre','_',...
            '10','percent/structs_p',num2str(pat),'pre.mat']; 
    
    currpat = load(filename_preStruct);
    n_tiles = 1;
    
    % Set up cell to save simulationr results
    SimOut = cell(n_stoch,n_tiles);
    
    for j = 1:size(param_array,1) % Parameter name loop
        
        pars_curr = param_array{j,2}; % take out params from param_array
        pars_name_curr = param_array{j,1}; % take out param name from param_array 
        disp(['Patient', num2str(pat),' par ',pars_name_curr]) % Display progress
        
        for ind_par_val = 1:size(pars_curr,2) % Parameter value loop
            
            % Get the current value for the parameter in question
            pars_curr_val = pars_curr(ind_par_val);
            
            % Set up data storage
            L_end = zeros(n_stoch,n_tiles,50); % L will store the RDF results
            L_lastobs = zeros(n_stoch,n_tiles,50);
            CD8_total_time = zeros(n_stoch,n_tiles,endTime); % Store time course of CD8 for all stoch runs for all tiles
            TUM_total_time = zeros(n_stoch,n_tiles,endTime); % Store time course of TUM for all stoch runs for all tiles
            CD8_total_time_tilemean = zeros(n_tiles,endTime); % Create matrices for mean of stoch runs for all tiles
            TUM_total_time_tilemean = zeros(n_tiles,endTime); % Create matrices for mean of stoch runs for all tiles
            L_end_tilemean = zeros(n_tiles,50); % for mean of RDF for each tile
            L_lastobs_tilemean = zeros(n_tiles,50); % for mean of RDF for each tile
            
            tile = 1;
            
                % Get the starting tile
                MySystem = currpat.allSystems{tiles_chosen(p)};
                
                % Inside the loop to run the simulations, save the full
                % output (every 2 timesteps) in SimOut (rows: stoch
                % repeats, cols:tiles)
                parfor stoch = 1:n_stoch % stochastic repeats
                    tic
                    curr = runSimulation(MySystem,stoch,pars_curr_val,pars_name_curr,0,0,0,endTime,2);% Last entry is 2 to indicate this is a sensitivity analysis
                    SimOut{stoch,tile}=curr;
                    toc
                end % end stoch repeats
            
            
            %% Get all the data out of the cell array so we can take
            % averages for RDF and cell numbers over time
            
               for stoch = 1:n_stoch
                   % If the simulation finished prematurely (tumour
                   % killed), indicate this in the RDF results
                    if size(SimOut{stoch,tile},2)<endTime
                        L_end(stoch,tile,:)=nan;
                        L_lastobs(stoch,tile,:)=computeRDF(SimOut{stoch,tile}{end}.grid.Li, 50);
                        
                        % Get the number of CD8 cells and tumour cells at each timepoint
                        for t_pt = 1:size(SimOut{stoch,tile},2)
                            CD8_total_time(stoch,tile,t_pt) = length(SimOut{stoch,tile}{t_pt}.IM.IMcells);
                            TUM_total_time(stoch,tile,t_pt) = length(SimOut{stoch,tile}{t_pt}.TU.TUcells);
                        end
                        
                        % If the simulations did not run to the end time,
                        % pad with NANs
                        for t_pt = size(SimOut{stoch,tile},2)+1:endTime
                            CD8_total_time(stoch,tile,t_pt) = nan;
                            TUM_total_time(stoch,tile,t_pt) = nan;
                        end
                    else
                        L_end(stoch,tile,:) = computeRDF(SimOut{stoch,tile}{endTime}.grid.Li, 50);
                        L_lastobs(stoch,tile,:) = L_end(stoch,tile,:);
                        
                        % Store number of CD8 and tumour cells at each time
                        % point
                        for t_pt = 1:endTime
                            CD8_total_time(stoch,tile,t_pt) = length(SimOut{stoch,tile}{t_pt}.IM.IMcells);
                            TUM_total_time(stoch,tile,t_pt) = length(SimOut{stoch,tile}{t_pt}.TU.TUcells);
                        end
                    end % if finished early loop
               end %stochastic loop
               
                %Put results into the matrices
                CD8_total_time_tilemean(tile,:) = nanmean(CD8_total_time(:,tile,:),1);
                TUM_total_time_tilemean(tile,:) = nanmean(TUM_total_time(:,tile,:),1);
                L_end_tilemean(tile,:) = nanmean(L_end(:,tile,:),1);
                L_lastobs_tilemean(tile,:)= nanmean(L_lastobs(:,tile,:),1);
            
                
            % Save all data for current paramter name, parameter value, current patient
            SA_allData{j,ind_par_val,p} = {pars_name_curr,pars_curr_val,pat,CD8_total_time_tilemean,TUM_total_time_tilemean,L_end_tilemean,L_lastobs_tilemean,CD8_total_time,TUM_total_time,L_end,L_lastobs};
        end % end param value
    end %end for param name
end % end for pats
mkdir(directory)
save(filenameout,'SA_allData')

