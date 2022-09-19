function run_make_exptFiles(animal, session, target_groups, saveDir, plane_spacing, this_dendrite,...
    num_blocks, iti, trials_per_block)

%% Dustin - 2021-01-22 (last update) 

num_trials_stim = trials_per_block; %photostim trials PER BLOCK, IMPORTANT: target_groups.num_trials has to add up to this value
num_trials_total = trials_per_block; %photostim trials plus blank trials PER BLOCK (num_trials_total-num_trials_stim = number of blank trials)

add_blank_flag = 0;
num_blocks = num_blocks;
stim_freq = 1/iti; % 1/isi in sec % was 1/1.23
jitter = 0;

% select if visual stimuli present - if yes, triggers need to be adapted
visual_stimuli_flag = 0;
offset_photoToVisual = 500-20; %how much to offset the photostim from visual stim triggers: note, the initial delay ...
% in MarkPoints execution has to be subtracted from the desired delay value (currently 20 ms) 
% if visual_stimuli_flag is 0, make offset_photoToVisual = NaN or 0;

% define number of targets and trials per pattern
num_trialTypes = size(target_groups,2);
trials_per_trialType = 1;

% define trial order
final_order = randsample(1:num_trialTypes,num_trialTypes);
final_order = reshape(final_order,[num_trials_stim,num_blocks])

% define individual pattern parameters 
num_unique_target_patterns = num_trialTypes;
cells_per_pattern = arrayfun(@(x) numel(x.these),target_groups);

%% MAKE EXPT FILES 
SLMtargets2 = {};SLMtargets2{1} = {};SaveNames2 = {};SaveNames2{1} = {};
% read yaml settings file to get some parameters
% just for DH office PC, run this: cd('C:\Users\Dustin\Dropbox\Bruker3')
yaml = ReadYaml('settings.yml');
ZeroOrderSLMCoordinates = [256 256];
SLMDimensions = [512, 512];
ZeroOrderSizePixels = yaml.ZeroOrderBlockSize_PX;

%% PREMAKE PHASE MASKS
allGP2 = NaN(num_unique_target_patterns,2);

if ~exist([saveDir filesep '999'], 'dir')
    mkdir([saveDir filesep '999'])
end
galvo_power_scale_factors = NaN(num_unique_target_patterns,1);

% just for DH office PC, run this: cd('C:\Users\Dustin\Dropbox\Bruker3\Naparm3\include')
for i = 1:num_unique_target_patterns
    
    disp(['computing galvo position for group ' num2str(i) ' of ' num2str(num_unique_target_patterns)])
    
    thisMask = target_groups(i).these_coordinatesXYZ; %target_groups(i).coordinates_selected_targets_XYZ;
    
    xx = thisMask(:,1);
    yy = thisMask(:,2);
    
    zz = NaN(numel(yy),1);
    for p = 1:numel(plane_spacing)
        this = plane_spacing(p);
        these = find(thisMask(:,3) == this);
        zz(these) = p;
    end
%     z_planes = unique(thisMask(:,3));
%     for p = 1:length(z_planes'
%         find(thisMask(:,3) == p)
%     end
%     
    thesePoints = [xx yy];  
    [OffsetPoints, GroupCentroid] = ZOBlockAvoider(thesePoints);
   
    allGP2(i,:) = GroupCentroid;
    SLMtargets2{1}{i} = [OffsetPoints(:,1) OffsetPoints(:,2) zz ones(length(zz),1)];
    
    SaveNames2{1}{i} = [...
        '999_'...
        num2str(i,'%04d')...
        '_' 'pattern' ...
        '_' num2str(length(zz)) 'Targets' ...
        '_X' num2str(GroupCentroid(1),'%03d') ...
        '_Y' num2str(GroupCentroid(2),'%03d') ...
        '.tif'];
    
    % new DH 2020-06-18: added weighting by galvo position
    galvo_power_scale_factors(i) = galvoPositionPowerWeighting(allGP2(i,1), allGP2(i,2));   
end

% make phase masks 
tic
[PhaseMasks, TransformedSLMTargets] = SLMPhaseMaskMakerCUDA3D_v3(...
    'Points', SLMtargets2{1},...
    'all_Galvo_Positions', allGP2, ...
    'plane_spacing', plane_spacing, ...
    'Save', true,...
    'SaveDirectory', [saveDir filesep '999'],...
    'SaveName', SaveNames2{1},...
    'Do3DTransform', true); % was false
toc 

premadePhaseMasks = PhaseMasks;
premadeTransformedSLMTargets = TransformedSLMTargets;
disp(['--- pre-made ' num2str(i) ' masks'])
clear PhaseMasks TransformedSLMTargets thisMask xx yy zz zumzum thesePoints OffsetPoints GroupCentroid

% check if all phasemasks look reasonable 
all_avg = [];
for q = 1:size(premadePhaseMasks,1)
    all_avg = [all_avg mean2(premadePhaseMasks{q})];
end
all_avg(find(isnan(all_avg))) = 0;
figure, histogram(all_avg)
title('check if any at 0. if yes, mask generation failed')
    
    
%% MAKE GPL, XML, MASKS AND TRIGGERS IN CORRECT ORDER FOR EACH BLOCK
SLMtargets = {};
% allX = nan(num_blocks, num_trials_stim, max([cells_per_pattern]));
% allY = nan(num_blocks, num_trials_stim, max([cells_per_pattern]));
% allZ = nan(num_blocks, num_trials_stim, max([cells_per_pattern]));
% offsetX = nan(num_blocks, num_trials_stim, max([cells_per_pattern]));
% offsetY = nan(num_blocks, num_trials_stim, max([cells_per_pattern]));
galvoX = nan(num_blocks, num_trials_stim, 1);
galvoY = nan(num_blocks, num_trials_stim, 1);
SaveNames = {};

% read yaml settings file to get some parameters
yaml = ReadYaml('settings.yml');
ZeroOrderSLMCoordinates = [256 256];
SLMDimensions = [512, 512];
ZeroOrderSizePixels = yaml.ZeroOrderBlockSize_PX;

% new DH 2020-01-22: add correction for group clustering
% clustered groups with many points close to zero order are more efficient
% than groups with more dispersed targets 
% used targetDispersion_compensationBruker3 to make fit. Load this fir here
% and apply inverse function to correct for this 

clust_powerCorrection = ...
    load(['C:\Users\User\Dropbox\Bruker3\PowerUtilities\dispersionTotalPowerCompensation_calFile2_smallClust.mat'],...
    'dispersion_px','powers','powers_norm','fittedmodel','dispersion_powerCal_mask')
patterns_chronological_allBlocks = [];

target_groups(1).cluster_correction_factor = [];
[target_groups.cluster_correction_factor] = deal(NaN);


for j = 1:num_blocks
    
    disp(['--- Designing experiment block ' num2str(j)])
    SLMtargets{j} = {};
    SaveNames{j} = {};
    disp('    Making targets and save names')
    
    thisBlockOrder = final_order(:,j);
    
    % iterate over trials for this block and make 
    % 1. galvo position list 
    % 2. file names 
    % 3. a record of all target coordinates 
    % for each trial
    
    patterns_chronological = NaN(num_trials_stim,1);
    counter_patterns = 1;
    
    for i = 1:num_trials_stim
        this_group = thisBlockOrder(i);
        
        thisMask = target_groups(this_group).these_coordinatesXYZ;
        
        xx = thisMask(:,1);
        yy = thisMask(:,2);
        zumzum = thisMask(:,3);
        
        thesePoints = [xx yy];  % note Y,X order
        [OffsetPoints, GroupCentroid] = ZOBlockAvoider(thesePoints);
        
        %         allX(j,i,1:size(thesePoints,1)) = thesePoints(:,1);
        %         allY(j,i,1:size(thesePoints,1)) = thesePoints(:,2);
        %         allZ(j,i,1:size(thesePoints,1)) = zumzum;
        %         offsetX(j,i,1:size(thesePoints,1)) = OffsetPoints(:,1);
        %         offsetY(j,i,1:size(thesePoints,1)) = OffsetPoints(:,2);
        galvoX(j,counter_patterns) = GroupCentroid(1);
        galvoY(j,counter_patterns) = GroupCentroid(2);
        
        SLMtargets{j}{counter_patterns} = [OffsetPoints(:,1) ...
            OffsetPoints(:,2) zumzum ones(length(xx),1)];  % x,y,z,I
        
        SaveNames{j}{counter_patterns} = [...
            num2str(j,'%02d') '_'...
            num2str(counter_patterns,'%04d')...
            '_' 'group' num2str(thisBlockOrder(i))...
            '_' num2str(length(xx)) 'Targets' ...
            '_X' num2str(galvoX(j,counter_patterns),'%03d') ...
            '_Y' num2str(galvoY(j,counter_patterns),'%03d') ...
            '.tif'];
        
        % compute degree of clustering (avg. offset from zero order)
        mean_disp = mean(pdist2(thesePoints, GroupCentroid));
        target_groups(this_group).mean_disp = mean_disp;
        
        cal_mask_disp = clust_powerCorrection.dispersion_powerCal_mask;
        this_mask_disp = mean_disp;
        target_groups(this_group).ref_disp = cal_mask_disp;
        
        correction_factor = ...
            clust_powerCorrection.fittedmodel(cal_mask_disp)/...
            clust_powerCorrection.fittedmodel(this_mask_disp);
        
        target_groups(this_group).cluster_correction_factor = correction_factor;

        patterns_chronological(counter_patterns) = this_group;
        counter_patterns=counter_patterns+1;
    end
   
    %% make GPL
    disp('    Making GPL')
    SpiralRevolutions = 3;
    SpiralDiameterUm  = 16;
    IsSpiral          = 'True';
    SaveName          = [saveDir filesep 'targetPerm_B' num2str(j)];
    %     MarkPoints_GPLMaker(galvoX(j,:), galvoY(j,:), IsSpiral, SpiralDiameterUm, SpiralRevolutions, SaveName);
    % on DH office PC run: cd('C:\Users\Dustin\Dropbox\Bruker3\MarkPoints')
    MarkPoints_GPLMaker_2(galvoX(j,:), galvoY(j,:), IsSpiral, ...
        SpiralDiameterUm, SpiralRevolutions, SaveName);
    % cd(saveDir)
    
    %% make XML
    disp('    Making XML')
    SaveName        = [saveDir filesep 'targetPerm_B' num2str(j)];
    NumGroups       = num_trials_stim; % multiply trial types by patterns per shot
    SequenceRepetitions = 1;
    NumRows         = NumGroups;
    ShotsPerPattern = 10;
    InterPatternShotInterval = 23;
    
    % define laser powers based on position and cell number in any given
    % stimulation group
    powers_per_group_MW = 18 * cells_per_pattern;%7.5
    LaserPowerMW      =  powers_per_group_MW(patterns_chronological) .* galvo_power_scale_factors(patterns_chronological)';
    
    % add correction for clustering 
    cluster_corr_factors = [target_groups(:).cluster_correction_factor];
    LaserPowerMW = LaserPowerMW .* cluster_corr_factors(patterns_chronological);
    LaserPowerPV    = round(mw2pv(LaserPowerMW));
    
%     % added DH 2021-02-11: correct first pulse power such that it's only
%     % 80% of the second pulse power to correct for heating related eff.
%     % drop-off 
%     heating_dropoff_vector = repmat([0.8 1], [1 numel(LaserPowerPV)/2]);
%     LaserPowerPV=LaserPowerPV.*heating_dropoff_vector;
    
    if any(LaserPowerPV > 1000)
        LaserPowerPV(find(LaserPowerPV > 1000)) = 1000;
        warning('max. stim laser power exceeded, set to 1000')
    end 
    
    InitialDelay    = 20;
    SpiralDuration  = 25;  % was 25 but with paired pulse stim you need an extra 20 ms for SLM update. trying to compensate for that by making spiral 1 ms shorter
    NumberOfTrials  = 1;
    AddDummy         = true;
    IterationDelay = 0;
    InterPointDelay = 24; % 0 = 100% duty cycle, 25 = 50% duty cycle 
    
    %max([0, InterPatternShotInterval-SpiralDuration]); 
    % note DH: ShotsPerPattern*(InterPointDelay+SpiralDuration) = total stim duration; add to this: InitialDelay
    
    % get parameters stored in settings file
    yaml = ReadYaml('settings.yml');
    VoltageOutputCategoryName = yaml.VoltageOutputCategoryName;
    VoltageOutputExperimentName = yaml.VoltageOutputExperimentName;
    LaserName = yaml.LaserName;
    TrigLine = yaml.TriggerLine;
        
    TriggerFreq                 = repmat({'First Repetition'}, NumRows, 1);
    TriggerSelect               = repmat({TrigLine}, NumRows, 1);
    AsyncSyncFrequency          = repmat({'FirstRepetition'}, NumRows, 1);
    VoltageOutputCategoryName   = repmat({VoltageOutputCategoryName}, NumRows, 1);
    VoltageOutputExperimentName = repmat({VoltageOutputExperimentName}, NumRows, 1);
    
    Indices = repmat(1:NumGroups,1,SequenceRepetitions)';
    PointNums = repmat(1:NumGroups,1,SequenceRepetitions)';
    Points = cell(NumRows,1);
    for p = 1:NumRows
        Points{p} = ['Point ' num2str(PointNums(p))];
    end
    
    Name = [...
        num2str(NumGroups) 'Patterns_' ...
        num2str(SequenceRepetitions) 'Repeats_x'...
        num2str(ShotsPerPattern) 'ShotPerPattern'...
        ];
    
    MarkPoints_XMLMaker(...
        'SaveName', SaveName, ...
        'ExptCat', 'NAPARM', ...
        'ExptName', Name, ...
        'NumRows', NumRows, ...
        'AddDummy', AddDummy, ...
        'UncagingLaser', LaserName,...
        'UncagingLaserPower', LaserPowerPV, ...
        'InternalIterations',SequenceRepetitions, ...
        'Repetitions', ShotsPerPattern, ...
        'InitialDelay', InitialDelay, ...
        'Duration', SpiralDuration, ...
        'InterPointDelay', InterPointDelay, ...
        'SpiralRevolutions', SpiralRevolutions, ...
        'TriggerFrequency', TriggerFreq, ...
        'TriggerSelection', TriggerSelect, ...
        'AsyncSyncFrequency', AsyncSyncFrequency, ...
        'VoltageOutputCategoryName', VoltageOutputCategoryName, ...
        'VoltageOutputExperimentName', VoltageOutputExperimentName, ...
        'Indices', Indices, ...
        'Points', Points, ...
        'Iterations', NumberOfTrials, ...
        'IterationDelay',IterationDelay ...
        );
    
    %% make triggers
    if ~add_blank_flag
        
        trigger_saveDir = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '\onlineMapping'];
        if ~exist(trigger_saveDir)
            mkdir(trigger_saveDir)
        end
        
        disp('    Making triggers')
        toblank = [];%find(all(isnan(selected_cells(theseIndices,:)),2)); % Lloyd used this to make stimuli blank that do not contain any more cells that are not in timeout
        AllTriggers = TriggerBuilder(...
            'name_prefix',        '',...
            'total_num_triggers', num_trials_stim,...
            'trigger_every_ms',   1000/stim_freq,...
            'trigger_dur_ms',     50,...
            'train_num_reps',     1,...
            'train_rep_every_ms', 0,...
            'shift_by_ms',        1,...
            'tail_to_add_ms',     0,...
            'jitter_ms',          jitter,...
            'trigger_amp_v',      5,...
            'sample_rate_hz',     1000,...
            'to_blank',           toblank,...
            'to_blank_train',     [],...
            'plot_result',        false,...
            'save_result',        false...
            );
        
        fid = fopen([trigger_saveDir filesep 'allPSTrig' num2str(j) '.dat'], 'w', 'l');
        fwrite(fid, AllTriggers, 'double');
        fclose(fid);
        
    elseif add_blank_flag
        % make all trial triggers - doesn't actually trigger anything but 
        disp('    Making triggers')
        toblank = [];%find(all(isnan(selected_cells(theseIndices,:)),2)); % Lloyd used this to make stimuli blank that do not contain any more cells that are not in timeout
        AllTriggers = TriggerBuilder(...
            'name_prefix',        '',...
            'total_num_triggers', num_trials_total,...
            'trigger_every_ms',   1000/stim_freq,...
            'trigger_dur_ms',     50,...
            'train_num_reps',     1,...
            'train_rep_every_ms', 0,...
            'shift_by_ms',        1,...
            'tail_to_add_ms',     0,...
            'jitter_ms',          jitter,...
            'trigger_amp_v',      5,...
            'sample_rate_hz',     1000,...
            'to_blank',           toblank,...
            'to_blank_train',     [],...
            'plot_result',        false,...
            'save_result',        false...
            );
        
%         trigger_saveDir = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\onlineMapping'];
        trigger_saveDir = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '\onlineMapping'];

        if ~exist(trigger_saveDir)
            mkdir(trigger_saveDir)
        end
        
        fid = fopen([trigger_saveDir filesep '2_allTrials_Triggers_B' num2str(j) '.dat'], 'w', 'l');
        fwrite(fid, AllTriggers, 'double');
        fclose(fid);
        
        % delete some triggers to generate blank stimuli 
        % ensure that
        % 1. first trial is not blank stimulus
        % 2. no two consecutive trials are blank (SLM might time out...)
        success = 0;
        while ~success
            blank_baselineStimuli = randsample(2:num_trials_total,num_trials_total-num_trials_stim);          
            diffs = diff(sort(blank_baselineStimuli));
            M = movsum(diffs,2);
            
            if any(M==2)
                success = 0;
            else 
                success = 1;
            end
        end
        
        blanks_all_blocks(:,j) = blank_baselineStimuli; % each column is a block
        
        AllTriggers_photo = zeros(length(AllTriggers),1);
        
        allTrig = thresholdDetect(AllTriggers,'above',4);
        if ~(length(allTrig) == num_trials_total)
            warning('error making triggers')
        end

        if ~visual_stimuli_flag
            for i = 1:length(allTrig)
                if ~ismember(i, blank_baselineStimuli)
                    thisI = allTrig(i);
                    thatI = thisI; %*randsample(0:7,1);
                    thoseI = [thatI:1:(thatI+49)];
                    AllTriggers_photo(thoseI) = 5;
                end
            end
        elseif visual_stimuli_flag % if visual stimuli present, offset visual and photostim by specific time delay
            for i = 1:length(allTrig)
                if ~ismember(i, blank_baselineStimuli)
                    thisI = allTrig(i);
                    thatI = thisI +  offset_photoToVisual/(1000/1000); %*randsample(0:7,1);
                    thoseI = [thatI:1:(thatI+49)];
                    AllTriggers_photo(thoseI) = 5;
                end
            end            
        end
        
        fid = fopen([trigger_saveDir filesep '1_toPV_Triggers_B' num2str(j) '.dat'], 'w', 'l');
        fwrite(fid, AllTriggers_photo, 'double');
        fclose(fid);

    end

    %% find (and save) phase masks from preMade batch
    
    disp('    Making phase masks')
    if ~exist([saveDir filesep num2str(j) filesep 'PhaseMasks'], 'dir')
        mkdir([saveDir filesep num2str(j) filesep 'PhaseMasks' ])
    end
    
    for mask = 1:length(patterns_chronological)
        
        thisPreMadeGroup = patterns_chronological(mask);
        thisSaveName = SaveNames{j}{mask};
        
        if isempty(strfind(thisSaveName,'.tiff')) && ~isempty(strfind(thisSaveName,'.tif'))
            imwrite(premadePhaseMasks{thisPreMadeGroup}, [saveDir filesep ...
                num2str(j) filesep 'PhaseMasks' filesep strrep(thisSaveName, '.tif', ['_CUDAphase.' 'tiff'])]);
        elseif ~isempty(strfind(thisSaveName,'.tiff'))
            disp('?')
            imwrite(premadePhaseMasks{thisPreMadeGroup}, [saveDir filesep ...
                num2str(j) filesep 'PhaseMasks' filesep strrep(thisSaveName, '.tiff', ['_CUDAphase.''tiff'])]);
        else
            sprintf('--- Warning: unknown image format')
        end 
    end

    patterns_chronological_allBlocks{j} = patterns_chronological;
end

powers_per_group_MW_scaled = galvo_power_scale_factors'.*powers_per_group_MW;
NumGroups = num_trialTypes;

cd(saveDir)
% save([saveDir filesep 'workspace.mat'])
save([saveDir filesep 'NAPARM_onlineSpineMapping_Points.mat'],...
    'target_groups')
if add_blank_flag
    save([saveDir filesep 'NAPARM_onlineSpineMapping_parameters.mat'],...
        'animal','session','num_trials_stim','num_trials_total','num_blocks','stim_freq','InitialDelay','InterPointDelay',...
        'SpiralDuration','AddDummy','IterationDelay','powers_per_group_MW','NumGroups',...
        'NumRows', 'ShotsPerPattern','SpiralRevolutions','final_order','yaml','allGP2','premadePhaseMasks',...
        'galvo_power_scale_factors','blanks_all_blocks','offset_photoToVisual','target_groups',...
        'patterns_chronological_allBlocks','num_trialTypes','trials_per_trialType','plane_spacing') 
else
    save([saveDir filesep 'NAPARM_onlineSpineMapping_parameters.mat'],...
        'animal','session','num_trials_stim','num_trials_total','num_blocks','stim_freq','InitialDelay','InterPointDelay',...
        'SpiralDuration','AddDummy','IterationDelay','powers_per_group_MW','NumGroups',...
        'NumRows', 'ShotsPerPattern','SpiralRevolutions','final_order','yaml','allGP2','premadePhaseMasks',...
        'galvo_power_scale_factors','offset_photoToVisual','target_groups',...
        'patterns_chronological_allBlocks','num_trialTypes','trials_per_trialType','plane_spacing') 
end

% % plot all groups 
% all_targets_figure = figure%('Position', [500 500 1000 1000]), 
% all_coords = vertcat(target_groups.these_coordinatesXYZ);
% scatter3(all_coords(:,1), all_coords(:,2), all_coords(:,3), 50, 'r', 'x'), hold on
% scatter3(this_soma_location(:,1), this_soma_location(:,2), mean(all_coords(:,3)), 75, 'b', '^')
% xlim([1 512]), ylim([1 512])
% set(gca,'YDir','reverse')
% legend('all targets','soma')
% 
% saveas(all_targets_figure,[saveDir filesep 'allTargetsFinal.fig']);


      
end