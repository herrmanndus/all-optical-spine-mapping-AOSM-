%% test sync paq for integrity

% last edit DH 2021-03-10: add option to read visual stim information


function sync_this_dendrite = extractPaq_dendri(animal, session, this_dendrite, numPlanes, nFrames, sync_saveDir, saveDir)

frameChannel = 'res_FrameStart';photostimChannel = 'toPV';photostimChannel_feeback = 'fromPV';
SLMchannel = 'fromSLM';allTrialsChannel = 'allTrials';photodiodeChannel = 'photoDiode';
% set dir
cd(sync_saveDir)

outputDir = saveDir;

% init output structure
sync_this_dendrite = {};

dendrite_id = this_dendrite;
sync_this_dendrite.dendrite = dendrite_id;

sync_fov = dendrite_id;
disp(['- dendrite ' num2str(dendrite_id) ' in FOV ' num2str(sync_fov)])

paq_files_thisD = dir([sync_saveDir '\' 'D' num2str(sync_fov) '*.paq']);
paq_files_thisD=natsortfiles({paq_files_thisD.name});

sync_this_dendrite.paq_files = paq_files_thisD;

photostim_files_dir = [saveDir];
stim_parameters_config = load([photostim_files_dir '\NAPARM_onlineSpineMapping_parameters.mat']);
% delete phase masks from the strucutre - slowing things down because large
% variable (1920,1180,420!) 
stim_parameters_config.premadePhaseMasks = [];

sync_this_dendrite.stim_parameters_config = stim_parameters_config;
sync_this_dendrite.animal = animal;
sync_this_dendrite.session = session;
sync_this_dendrite.sync_fov = sync_fov;

for block = 1:length(paq_files_thisD)
    
    [paqdata1, paqchannels] = paq2lab([sync_saveDir '\' paq_files_thisD{1,block}]);
    
    frameTimes1 = thresholdDetect(paqdata1(:,ismember(paqchannels, frameChannel)), ...
        'above', 1);
    
    nframes1 = nFrames;
    frameTimes1 = frameTimes1(1:nframes1);
    isiTooBig = find(diff(frameTimes1) > 10*median(diff(frameTimes1)));
    frameTimes1((isiTooBig+1):1:length(frameTimes1)) = [];
    frameTimes1Planes = reshape(frameTimes1, numPlanes, nframes1/numPlanes);
    
    % to PV
    photostimTimes=[];
    framesWithPhotostim =[];
    for i = 1:numPlanes
        edges = [frameTimes1Planes(i,1:end) +Inf];
        photostimTimes(i,:) = thresholdDetect(paqdata1(:,ismember(paqchannels, photostimChannel)), 'above', 1);
        framesWithPhotostim(i,:) = discretize(photostimTimes(i,:), edges);
    end
    
    % feedback from PV
    photostimTimes_fb=[];
    framesWithPhotostim_fb =[];
    for i = 1:numPlanes
        %         edges = [-Inf, mean([frameTimes1Planes(i,2:end)' frameTimes1Planes(i,1:end-1)'],2)', +Inf];
        edges = [frameTimes1Planes(i,1:end) +Inf];
        photostimTimes_fb(i,:) = thresholdDetect(paqdata1(:,ismember(paqchannels, photostimChannel_feeback)), 'above', 1);
        framesWithPhotostim_fb(i,:) = discretize(photostimTimes_fb(i,:), edges);
    end
    
    %         short_isi = diff(photostimTimes_fb); short_isi=short_isi(find(short_isi < mean(short_isi)));
    %         figure, histogram(short_isi)
    
    
    sync_this_dendrite.sync_blocks{block}.name = ...
        paq_files_thisD{1,block};
    sync_this_dendrite.sync_blocks{block}.block = ...
        block;
    sync_this_dendrite.sync_blocks{block}.framesWithPhotostim = ...
        framesWithPhotostim;
    sync_this_dendrite.sync_blocks{block}.frameTimes1Planes = ...
        frameTimes1Planes;
    sync_this_dendrite.sync_blocks{block}.nframes1 = ...
        nframes1;
    sync_this_dendrite.sync_blocks{block}.photostimTimes = ...
        photostimTimes;
    
    if size(framesWithPhotostim_fb,2) ==  size(framesWithPhotostim,2) & ...
            sum(sum(framesWithPhotostim-framesWithPhotostim_fb,[])) < 10 & ...
            size(framesWithPhotostim_fb,2) == sync_this_dendrite.stim_parameters_config.num_trials_stim
        disp(['-- feedback from PV ok'])
        disp(['-- ' num2str(size(framesWithPhotostim_fb,2)) ' stim/block'])
%         disp(['-- ' num2str(size(stim_parameters_points.target_groups,2)) ' stim groups'])
    else
        warning(['feedback from PV possibly NOT ok, please check manually'])
    end
    
    % also check SLM
    SLM_times=[];
    SLM_frames =[];
    for i = 1:numPlanes
        %         edges = [-Inf, mean([frameTimes1Planes(i,2:end)' frameTimes1Planes(i,1:end-1)'],2)', +Inf];
        edges = [frameTimes1Planes(i,1:end) +Inf];
        SLM_times(i,:) = thresholdDetect(paqdata1(:,ismember(paqchannels, SLMchannel)), 'above', 1);
        SLM_frames(i,:) = discretize(SLM_times(i,:), edges);
    end
    if size(framesWithPhotostim_fb,2) == length(SLM_frames(i,:))
        disp(['-- SLM updates ok'])
    elseif size(framesWithPhotostim_fb,2) < length(SLM_frames(i,:))
        warning('extra SLM triggers, go in and check!')
        pause
    elseif size(framesWithPhotostim_fb,2) > length(SLM_frames(i,:))
        warning('SLM triggers missing, go in and check!')
        pause
    elseif length(SLM_frames(i,:)) == 0
        warning('no SLM updates, check if intentional!')
        pause
    else
        disp(['-- SLM update error, please check manually'])
        pause
    end
    
    % get 'allTrials' times. These signal which trials were photostim or 'blank' trials
    % feedback from PV
    allTrials_times =[];
    allTrials_times_original =[];
    allTrials_frames =[];
    allTrials_frames_original =[];
    for i = 1:numPlanes
        %         edges = [-Inf, mean([frameTimes1Planes(i,2:end)' frameTimes1Planes(i,1:end-1)'],2)', +Inf];
        edges = [frameTimes1Planes(i,1:end) +Inf];
        allTrials_times(i,:) = thresholdDetect(paqdata1(:,ismember(paqchannels, allTrialsChannel)), 'above', 1);
        allTrials_times_original=allTrials_times;
        
        allTrials_frames(i,:) = discretize(allTrials_times(i,:), edges);
        allTrials_frames_original(i,:) = discretize(allTrials_times_original(i,:), edges);
        
    end  
           
    [frames_ps,trials_ps] = intersect(allTrials_frames(i,:),framesWithPhotostim(i,:));
    [frames_blank,trials_blank] = setxor(allTrials_frames(i,:),framesWithPhotostim(i,:));
    
    if sync_this_dendrite.stim_parameters_config.num_trials_total == numel(allTrials_frames) & ...
            sync_this_dendrite.stim_parameters_config.num_trials_total-sync_this_dendrite.stim_parameters_config.num_trials_stim == numel(trials_blank)
        disp('-- blanks detected')
    else
        frames_ps = [];trials_ps = [];frames_blank = [];trials_blank = [];
        for t=1:numel(allTrials_frames)
            this = allTrials_frames(t);these = [this-1:this+1];
            if sum(ismember(these,framesWithPhotostim)) == 1
                frames_ps = [frames_ps this];
                trials_ps = [trials_ps; t];
            else
                frames_blank = [frames_blank this];
                trials_blank = [trials_blank; t];
            end
        end
        if sync_this_dendrite.stim_parameters_config.num_trials_total == numel(allTrials_times) & ...
                sync_this_dendrite.stim_parameters_config.num_trials_total-sync_this_dendrite.stim_parameters_config.num_trials_stim == numel(trials_blank)
            disp('-- blanks detected')
        else
            warning('-- wrong number of blanks, please check manually')
            pause
        end
    end
    
    sync_this_dendrite.sync_blocks{block}.frames_ps = frames_ps;
    sync_this_dendrite.sync_blocks{block}.trials_ps = trials_ps;
    sync_this_dendrite.sync_blocks{block}.frames_blank = frames_blank;
    sync_this_dendrite.sync_blocks{block}.trials_blank = trials_blank;
    sync_this_dendrite.sync_blocks{block}.allTrials_frames_original = allTrials_frames_original;
    
%     % check for correct ps/blank order
%     if ~numel(intersect(trials_blank,sync_this_dendrite.stim_parameters_config.blanks_all_blocks(:,block))) == numel(trials_blank)
%         warning('blank order wrong, check manually')
%         pause
%     else
%         disp(['-- ' num2str(numel(frames_blank)) ' blanks/block'])
%     end
    
    disp(['--- block ' num2str(block) ' of ' num2str(length(paq_files_thisD)) ' done']) 
end

save([saveDir filesep 'syncAllDendrites_M' num2str(animal) '_S' num2str(session)   '.mat'], 'sync_this_dendrite', 'animal', 'session', '-v7.3')










% for ... extract all tuning relevant vaiables and store in second row of output strucutre
% % find all tuning/visual stim sessions without photostim
% tuningFiles = dir([directory 'visual*\tuning*']);
% allFilesSortedVis = natsortfiles({tuningFiles.name});
% end


