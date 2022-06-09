%% Oliver Grimm and Lucy Hutchinson 2021
%
function [] = runsims_parsets(pars_all, par_names, pats_chosen, endtimes, sim_endtimes,...
    tilenumbers, tile_limit,  n_stoch, save_path, task_name,savetimecourse,saveselectedgridstack)%

n_pars = length(pars_all);

%% 2. Run simulations for chosen patients for every parameter sample and save the full outputs

for p = 1:length(pats_chosen)
    pat = pats_chosen(p);
    obs_endTime = endtimes(p); % Observed timepoint on treatment biopsy
    sim_endTime = sim_endtimes(p); % Time to which we wish to simulate (e.g. endTime+10)
    n_tiles_curr = tilenumbers(p);
    
    % For the patients with more than tile_limit tiles, we want to
    % subsample tile_limit
    % tiles at random to improve computation time and avoid massive
    % overrepresentation of some patients
    if n_tiles_curr > tile_limit && ~isempty(tile_limit)
        tile_inds = randsample(n_tiles_curr,tile_limit,'false');
        tile_inds = sort(tile_inds);
        n_tiles_sim=tile_limit;
    else
        tile_inds = 1:n_tiles_curr;
        n_tiles_sim = n_tiles_curr;
    end    
    
    % Set up a vector of the instance to save [pat, tile, parset, stoch] if
    % applicable
    instance_to_save = [];
    if ~isempty(saveselectedgridstack)
            
        if ismember(pat,saveselectedgridstack(:,1))
            % Extract the row of details of what to save
            instance_to_save = saveselectedgridstack(saveselectedgridstack(:,1)==pat,:);
            
            % In case we try to save a tile that doesn't exist, save the
            % last one (save for parsets and stoch runs)
            if instance_to_save(2)>n_tiles_sim
                instance_to_save(2)=n_tiles_sim;
            end
            if instance_to_save(3)>n_pars
                instance_to_save(3)=n_pars;
            end
            if instance_to_save(4)>n_stoch
                instance_to_save(4)=n_stoch;
            end
       
        end
        
    end
    
    
    % read in pre- structs
   filename_preStruct = ['final_tiles/PatientData/p',num2str(pat),'/','pre','_',...
        '10','percent/structs_p',num2str(pat),'pre.mat'];
    
    current_pat = load(filename_preStruct);
        
    % Set up storage
    nCD8_endtime_all={};
    nTUM_endtime_all={};
    flags_finished_all={};
    earlyEndtime_all={};
    CD8gridstack_endtime_all={};
    TUMgridstack_endtime_all={};
    t_currpars_currpat=[];
    CD8_timecourse_all=nan(n_tiles_sim,sim_endTime,n_stoch,size(pars_all,1));
    TUM_timecourse_all=nan(n_tiles_sim,sim_endTime,n_stoch,size(pars_all,1));
    simout_selected = {}; % for saving one chosen instance
   
    % Run simulations for chosen patient for every parameter sample and
    % every stochastic run and every tile and save the results at the last
    % observed time point
    parfor j = 1:n_pars % Parameter set loop
        
        pars_curr = pars_all(j,:);
        
        t_currpars = tic; % save the time taken for each par set
        
        % Set up the output storage
        SimOut = cell(n_stoch,n_tiles_curr);
        
        for tt = 1:n_tiles_sim
            tile = tile_inds(tt); % for each tile
            MySystem = current_pat.allSystems{tile};
            
            for stoch = 1:n_stoch
                disp(['Patient ', num2str(pat),' sample ',num2str(j),', tile ',num2str(tile),', nstoch ',num2str(stoch)])
                tic
                SimOut{stoch, tt} = runSimulation(MySystem,stoch,pars_curr,par_names,0,0,0,sim_endTime,1);
                toc
            end % end stoch repeats
            
        end % end tile loop
        
        % Put the results into cell arrays so all results for each patient
        % are together
        nCD8_endtime=zeros(n_tiles_sim,n_stoch);
        nTUM_endtime=zeros(n_tiles_sim,n_stoch);
        flags_finished = zeros(n_tiles_sim,n_stoch);
        earlyEndtime = zeros(n_tiles_sim,n_stoch);
        CD8gridstack_endtime = cell(n_tiles_sim,n_stoch);
        TUMgridstack_endtime = cell(n_tiles_sim,n_stoch);
        CD8_timecourse = nan(n_tiles_sim,sim_endTime,n_stoch);
        TUM_timecourse = nan(n_tiles_sim,sim_endTime,n_stoch);
        
        for tt = 1:n_tiles_sim
            for stoch = 1:n_stoch
                earlyEndtime(tt,stoch)=size(SimOut{stoch,tt},2);
                timepoint_curr = min(earlyEndtime(tt,stoch), obs_endTime);
                nCD8_endtime(tt,stoch) = sum(sum(SimOut{stoch,tt}{1,timepoint_curr}.grid.Li));
                nTUM_endtime(tt,stoch) = sum(sum(SimOut{stoch,tt}{1,timepoint_curr}.grid.Lt));
                CD8gridstack_endtime{tt,stoch}=logical(SimOut{stoch,tt}{1,timepoint_curr}.grid.Li);
                TUMgridstack_endtime{tt,stoch}=logical(SimOut{stoch,tt}{1,timepoint_curr}.grid.Lt);
                % If we want to save timecourses of cell numbers, this is
                % where we do it
                if savetimecourse
                    for timept = 1:length(SimOut{stoch, tt})
                        nCD8_curr = sum(sum(SimOut{stoch, tt}{timept}.grid.Li));
                        nTUM_curr = sum(sum(SimOut{stoch, tt}{timept}.grid.Lt));
                        CD8_timecourse(tt,timept,stoch)=nCD8_curr;
                        TUM_timecourse(tt,timept,stoch)=nTUM_curr;
                    end
                end
                
                % Extract the selected instance 
                
                if ~isempty(instance_to_save) 
                    if instance_to_save == [pat, tt, j, stoch]
                    simout_selected = SimOut{stoch,tt};    
                    end
                end  
                                
                if earlyEndtime(tt,stoch)<obs_endTime
                    flags_finished(tt,stoch) = 0;
                else
                    flags_finished(tt,stoch) = 1;
                end
                
            end
        end
        
        nCD8_endtime_all{1,j}=nCD8_endtime;
        nTUM_endtime_all{1,j}=nTUM_endtime;
        flags_finished_all{1,j}=flags_finished;
        earlyEndtime_all{1,j}=earlyEndtime;
        CD8gridstack_endtime_all{1,j}=CD8gridstack_endtime;
        TUMgridstack_endtime_all{1,j} = TUMgridstack_endtime;
        t_currpars_currpat(1,j)=toc(t_currpars);
        CD8_timecourse_all(:,:,:,j) = CD8_timecourse;
        TUM_timecourse_all(:,:,:,j) = TUM_timecourse;
        
        
    end % end pars
    % save here
    save([save_path, '/',task_name,'_pat',num2str(pat)],'tile_inds','nCD8_endtime_all','nTUM_endtime_all',...
        'flags_finished_all','earlyEndtime_all','CD8gridstack_endtime_all','TUMgridstack_endtime_all',...
        't_currpars_currpat','CD8_timecourse_all','TUM_timecourse_all', 'tile_inds','simout_selected','instance_to_save');
end % end for pats

