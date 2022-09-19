% AOSM grid maker script 
% Author: Dustin Herrmann, 2021-2022
% allows making semi-automated grid to stimulate LM locations while imaging
% dendrite in V1
% for making of synatpical connections using the all-optical method

%% DEFINE SESSION PARAMETERS

animal = 1095;
session = 3;
widefield_session = 1;

current_session_dir = ['E:\Data\M' num2str(animal) '\S' num2str(session)];
widefield_session_dir = ['Y:\M' num2str(animal) '\S' num2str(widefield_session)];

saveDir_allDendrites = ['E:\Data\M' num2str(animal) '\S' num2str(session) filesep 'targetPerm_onlineMapping']; 

if ~exist(saveDir_allDendrites)
    mkdir(saveDir_allDendrites)
end

%% LOAD COORDINATES AND SURFACE IMAGES
% target_candidates = load([current_session_dir '\export_points.mat']);

% load 2-photon FOV image of the braun surface with clearly visible blood
% vessel pattern
temp_surface_file = dir([current_session_dir filesep '\*\' 'surface-001_Cycle00001_Ch2_000001.ome.tif']);
current_surface = uint8(imread([temp_surface_file.folder filesep temp_surface_file.name])); clear temp_surface_file

% load widefield FOV image  
widefield_surface = load([widefield_session_dir filesep 'INTdualMonMaps.mat'],'FOVGiRegRot');


%% MAKE QUICK TRANSFORM 
% manually select matching coordinates based on blood vessels to compute a
% transform from 

[selectedMovingPoints,selectedFixedPoints] = cpselect(current_surface,widefield_surface.FOVGiRegRot,'Wait',true)
tform = fitgeotrans(selectedMovingPoints,selectedFixedPoints,'affine')

% plot to manually check if transform does what it should 
Jreg = imwarp(current_surface,tform,'OutputView',imref2d(size(widefield_surface.FOVGiRegRot)));
reg_test_figure = figure('Position', [500 500 1500 1500]), hold on 
imshowpair(widefield_surface.FOVGiRegRot,Jreg,'blend','Scaling','joint')
title(['M: ' num2str(animal), ' S: ' num2str(session) ', WF: ' num2str(widefield_session)])

% save figure and transform details 
saveas(reg_test_figure,[saveDir_allDendrites filesep 'surface_registration.fig']);
save([saveDir_allDendrites filesep 'thisSessionVariables_M' num2str(animal) '_S' num2str(session) '.mat'],...
    'animal','session','widefield_session','current_session_dir','widefield_session_dir',...
    'saveDir_allDendrites','current_surface','widefield_surface','selectedMovingPoints','selectedFixedPoints',...
    'tform','Jreg')



%% LOAD Z-STACK AND ALLOW MANUAL SELECTION OF SOMA LOCATION
% load Z-stack of 2-photon fov size. select position of the soma. then transform between 2-photon and widefield space 

% define dendrite ID 
this_dendrite =1;

% load stack file and select soma location
stack_file_dir = dir([current_session_dir filesep '*/*fullStack-001*_Ch2.tif']);
stack_file = uint8(TiffReaderAP([stack_file_dir.folder filesep stack_file_dir.name])); clear stack_file_dir

StackSlider(stack_file), disp('select appropriate slice, press enter in command window when ready')
pause

thisROI=imfreehand
this_soma_location = mean(getPosition(thisROI),1);close

% transform soma location and target candiates to WF space and get azi/elev preferences 
[soma_x, soma_y] = transformPointsForward(tform, this_soma_location(1), this_soma_location(2));
WF_soma_px = round([soma_x soma_y]); clear soma_x soma_y


%% INFER RETINOTOPIC LOCATIONS FROM WF MAP AND SPLIT TARGETS INTO NEAR AND FAR  

% load widefield retinotopic map (if required) and get soma azimuth and
% elevation locations 

LR_map = load([widefield_session_dir '\INTdualMonMaps.mat'],'aziMapfilt3','aziSpan','aziMapRaw');
UD_map = load([widefield_session_dir '\INTdualMonMaps.mat'],'elevMapfilt3','elevSpan','elevMapRaw');

LR_map.aziMapfilt5 = medfilt2(LR_map.aziMapRaw,[5 5]);UD_map.elevMapfilt5 = medfilt2(UD_map.elevMapRaw,[5 5]);

mid_line_phase = LR_map.aziSpan(1,1);lateral_phase =  LR_map.aziSpan(1,2);mid_line_deg = LR_map.aziSpan(2,1);lateral_deg =  LR_map.aziSpan(2,2);
up_phase = UD_map.elevSpan(1,1);down_phase = UD_map.elevSpan(1,2);up_deg = UD_map.elevSpan(2,1);down_deg = UD_map.elevSpan(2,2);

LR_map.degrees = interp1([mid_line_phase,lateral_phase],...
    [mid_line_deg,lateral_deg],LR_map.aziMapfilt5,'linear','extrap');
UD_map.degrees = interp1([up_phase,down_phase],...
    [up_deg,down_deg],UD_map.elevMapfilt5,'linear','extrap');


%% MAKE STIMULATION GRID

saveDir = [saveDir_allDendrites filesep 'D' num2str(this_dendrite)];

%%%%%%%
grid_spacing = 23; % in um

[images_allPlanes,grid_manual_selection] = ...
    makeTargetGrid_retEx3(current_session_dir,grid_spacing,saveDir,tform,this_soma_location,WF_soma_px,...
    LR_map, UD_map); % --> triangular grid 

%%%%%%%
retinotopic_range = [0 50];
plane_spacing = [-200 -150 -100 -50 0 50]; % number of planes has to match the number of c1v1 mean images 
set_ae_sep_flag = 0; % set to 1 if you want to apply separate azimuth and elevation ranges using azi_range & elev_range
azi_range = [0 25]; % if set_ae_sep_flag=1. use this as the azimuth range 
elev_range = [0 25]; % if set_ae_sep_flag=1. use this as the elevation range 
shift_soma_AE = [0 0]; % by how much to shift soma azimuth and elevation position to account for LM mapping issues 
[grid_coordinates_planesXYZP, grid_coordinates_planesWFxy_AE] = ...
    refineGridRet(images_allPlanes,grid_manual_selection,tform,LR_map,UD_map,WF_soma_px,retinotopic_range,plane_spacing,saveDir,this_soma_location, azi_range, elev_range, set_ae_sep_flag, shift_soma_AE);

%%%%%%%
n_targets = 25;
num_blocks = 3; % was 2 
trials_per_block = 420; % note: if you change this, you need to change 'iti' accordingly (see below) 
[target_groups]=makeStimGroupsGrid(grid_coordinates_planesXYZP,n_targets,saveDir, num_blocks, trials_per_block);
%%%%%%%

%%%%%%%
iti = 1.23; % in seconds 


%% Generate all files needed for the experiment 
run_make_exptFiles(animal, session, target_groups, saveDir, plane_spacing, this_dendrite, num_blocks, iti, trials_per_block)
% change power line 238

%% plot some stuff and make ret reference 

% get soma ret location 
[A,B] = meshgrid(WF_soma_px(2)-1:WF_soma_px(2)+1,WF_soma_px(1)-1:WF_soma_px(1)+1);soma_halo = [A(:) B(:)];
soma_ind = sub2ind(size(LR_map.degrees),soma_halo(:,1),soma_halo(:,2));
soma_AE = [median(LR_map.degrees(soma_ind)), median(UD_map.degrees(soma_ind))];

[target_x, target_y] = transformPointsForward(tform, grid_coordinates_planesXYZP(:,1), grid_coordinates_planesXYZP(:,2));
WF_target_px = round([target_x target_y]); clear target_x target_y

grid_coordinates_planesXYZP_WFXY = [grid_coordinates_planesXYZP WF_target_px];

% get all target candidate ret locations 
candidates_AE = NaN(size(WF_target_px,1),2);
for cand = 1:size(WF_target_px,1)
    this_px = WF_target_px(cand,:);
    [A,B] = meshgrid(this_px(2)-1:this_px(2)+1,this_px(1)-1:this_px(1)+1);target_halo = [A(:) B(:)];
    target_ind = sub2ind(size(LR_map.degrees),target_halo(:,1),target_halo(:,2));

    target_AE = [median(LR_map.degrees(target_ind)), median(UD_map.degrees(target_ind))];
    candidates_AE(cand,:) = target_AE;
end 
    
AE_scatter = figure('Position', [500 500 1000 1000]), scatter(soma_AE(1), soma_AE(2),'b','filled'), hold on, scatter(candidates_AE(:,1), candidates_AE(:,2),'r','filled')
xlabel('Azimuth'), ylabel('Elevation'), legend('soma','targets'), saveas(AE_scatter,[saveDir filesep 'AE_scatter.fig']);

grid_coordinates_planesXYZP_WFXY_AE = [grid_coordinates_planesXYZP_WFXY candidates_AE];

% plot WF map and save it 

WF_target_locations = figure, 
WF_target_locations.Position = [800 400 2000 800];

subplot(121),
    imagesc(LR_map.degrees), hold on,
    axis off
    scatter(WF_target_px(:,1), WF_target_px(:,2), 20, 'k', '*')
    scatter(WF_soma_px(1), WF_soma_px(2), 100, 'b', 'filled', '^')
    title('Azimuth - all targets')
subplot(122),
    imagesc(UD_map.degrees), hold on,
    axis off
    scatter(WF_target_px(:,1), WF_target_px(:,2), 20, 'k', '*')
    scatter(WF_soma_px(1), WF_soma_px(2), 100, 'b', 'filled', '^')
    title('Elevation - all targets')
saveas(WF_target_locations,[saveDir filesep 'WF_target_locations.fig']);


% plot to verify transforms 
all_candidates_plot = figure%('Position', [500 500 2200 1000]), 
subplot(121), imagesc(current_surface), colormap(gray), hold on, 
scatter(this_soma_location(1), this_soma_location(2),100,'w','filled')
scatter(grid_coordinates_planesXYZP(:,1), grid_coordinates_planesXYZP(:,2), 20, 'r', '*')

subplot(122), imagesc(widefield_surface.FOVGiRegRot), colormap(gray), hold on, 
scatter(WF_target_px(:,1), WF_target_px(:,2), 20, 'r', '*')
scatter(WF_soma_px(1), WF_soma_px(2),100,'w','filled'), 

saveas(all_candidates_plot,[saveDir filesep 'all_candidates_physLocation.fig']);

twoP_target_locations = figure,    
twoP_target_locations.Position = [800 400 2000 800];
subplot(121),
    imagesc(current_surface), hold on, 
    colormap(gray)
    axis off
    coords = vertcat(target_groups(:).these_coordinatesXYZ);
    scatter(grid_coordinates_planesXYZP(:,1), grid_coordinates_planesXYZP(:,2), 20, 'r', '*')
    scatter(this_soma_location(1), this_soma_location(2), 100, 'w', 'filled', '^')
    title('all target locations')
subplot(122),
    imagesc(current_surface), hold on, 
    colormap(gray)
    axis off
    random_selection = randsample(1:size(target_groups,2),5);
    cmp = jet(length(random_selection));
    for example = 1:length(random_selection)
        coords = target_groups(random_selection(example)).these_coordinatesXYZ;
        scatter(coords(:,1), coords(:,2), 20,cmp(example,:),  '*')
    end
    scatter(this_soma_location(1), this_soma_location(2), 100, 'w', 'filled', '^')
    title([num2str(length(random_selection)) ' example groups'])

% saveas(WF_target_locations,[saveDir filesep 'WF_target_locations.fig']);
saveas(twoP_target_locations,[saveDir filesep 'twoP_target_locations.fig']);


%% SET EXPT PARAMETERS AND MAKE TRIAL ORDER 
save([saveDir filesep 'targetPerm_ONLINE_details_M' num2str(animal) '_S' num2str(session) '_D' num2str(this_dendrite) '.mat'], ...
    'animal', 'session','widefield_session','current_session_dir','widefield_session_dir',...
    'saveDir','current_surface','widefield_surface',...
    'selectedMovingPoints','selectedFixedPoints','tform','this_dendrite','stack_file',...
    'this_soma_location','WF_soma_px','WF_target_px','LR_map','UD_map',...
    'soma_AE','candidates_AE','target_groups',...
    'grid_coordinates_planesXYZP_WFXY_AE','plane_spacing','grid_spacing','n_targets','retinotopic_range')

save([saveDir filesep 'workspace_M' num2str(animal) '_S' num2str(session) '_D' num2str(this_dendrite) '.mat'])


%% ONLINE REGISTRATION SETUP

% finalize FOV, take 500-1000 frames at 1024x1024 full resolution
% call it regref_D*, enable conversion to tiff!!

% run this: and select the regref_D folder, takes about 30 seconds
stack_mean_im(1)

% go to folder regref_D*-00*_stack and inspect registered image - should look like registered file of
% current FOV --> this is your registration reference frame

PrairieLink_RawDataStreamReg
% click "Do registration", select registered regref image 

% disable conversion again! 

%% NOW RUN ONLINE MAPPING EXPERIMENT - 3 blocks! 
% phase masks and galvo position lists are in this folder: targetPerm_onlineMapping

%% RUN SYNC EXTRACTION 
% sync_saveDir = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\onlineMapping'];
sync_saveDir = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '\onlineMapping'];

nFrames = 8000;
tic
sync_this_dendrite = extractPaq_dendri(animal, session, this_dendrite, 1, nFrames, sync_saveDir, saveDir)
toc
% takes ~30 seconds 


%% make target groups with online mapped targets 

saveDir_allDendrites2 = ['E:\Data\M' num2str(animal) '\S' num2str(session) filesep 'targetPerm_sigTargets']; 

if ~exist(saveDir_allDendrites2)
    mkdir(saveDir_allDendrites2)
end

saveDir2 = [saveDir_allDendrites2 filesep 'D' num2str(this_dendrite)]; 

if ~exist(saveDir2)
    mkdir(saveDir2)
end

% edit onlineStimResponseMapping
edit online_spine_detection


% compute target dispersion for power compensation
% compute dispersion
for mask = 1:size(enrichedTargets,2)
    these_coordinatesXYZ = enrichedTargets(mask).these_coordinatesXYZ;
    [~, GroupCentroid] = ZOBlockAvoider(these_coordinatesXYZ(:,1:2));
    this_disp = mean(pdist2(these_coordinatesXYZ(:,1:2),GroupCentroid));
    enrichedTargets(mask).this_disp=this_disp;
end


run_make_exptFiles_candidates(animal, session, enrichedTargets, this_soma_location, saveDir2, this_dendrite)





%%










%% new section: post hoc check if global faciliation 
% to read sync from sparse noise blocks run this 
sync_saveDir2 = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '\sigTargets_SN'];
nFrames2=4000;
sync_posthoc = extractPaq_dendrites_targetPermPOSTHOC_SN(animal, session, this_dendrite, 1, nFrames2, sync_saveDir2, saveDir2)

% to read sync from annulus blocks run this 
sync_saveDir2 = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '\sigTargets_ANN'];
sync_posthoc = extractPaq_dendrites_targetPermPOSTHOC_ANN(animal, session, this_dendrite, 1, nFrames, sync_saveDir2, saveDir2)

% for older sessions (e.g. darkness run this
sync_saveDir2 = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '\sigTargets'];
sync_posthoc = extractPaq_dendrites_targetPermPOSTHOC(animal, session, this_dendrite, 1, nFrames, sync_saveDir2, saveDir2)

% GLOBAL
allSP=readtable(['E:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) 'globalF.csv']);
thisF=table2array(allSP(:,[2]))';

figure, plot(thisF)%hold on; plot(thisF2); figure,plot(thisF-thisF2);thisF1=thisF;thisF=thisF-thisF2;

ps_times = [];
for block = 1:numel(sync_posthoc.sync_blocks)
    these = sync_posthoc.sync_blocks{block}.framesWithPhotostim;
    ps_times = [ps_times (these+(block-1)*8000)]
end

% ps_times = floor(ps_times/4);

all_trials = [];
for ps = 1:length(ps_times)
    this_f = ps_times(ps);
    globalF = thisF((this_f-8):(this_f+16));
    all_trials = [all_trials;globalF];
end
    
global_figure = figure, plot(all_trials'), hold on, plot(mean(all_trials),'Color','k','LineWidth',4)
yy = ylim;
line([9 9],[yy])
line([11 11],[yy])

saveas(global_figure,[saveDir_allDendrites2 filesep 'D' num2str(this_dendrite) '_global_posthocSTA.fig']);


% HANDPICKED SPINES 
allSP=readtable(['E:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '_posthoc_spines.csv']);
% allSP=readtable(['Y:\M1079\S5' '\D' num2str(this_dendrite) '_posthoc_spines.csv']);
this=table2array(allSP(:,[3:4:size(allSP,2)]));
numSpines=size(this,2);
save(['D' num2str(dendrite) 'numSpines'], 'numSpines');

spinesF=NaN(size(this,1),numSpines); 
for spine=1:numSpines
    thisTrace=this(:,spine);
    thisF0=running_percentile(thisTrace,1350,10);%thisF0=thisF0';
    thisdF=thisTrace-thisF0; thisdFF=thisdF./thisF0;
    spinesF(:,spine) = thisdFF;
    fprintf('spine %0d F obtained\n',spine)
end

spinesF = spinesF(:,handpicked);

for sp = 1:size(spinesF,2)
    
    thisF = spinesF(:,sp)';
    
    ps_times = [];
    for block = 1:numel(sync_posthoc.sync_blocks)
        these = sync_posthoc.sync_blocks{block}.framesWithPhotostim;
        ps_times = [ps_times (these+(block-1)*8000)]
    end
    %     ps_times = ps_times(1:55)
    %     ps_times = floor(ps_times/4);
    
    all_trials = [];
    for ps = 1:length(ps_times)
        this_f = ps_times(ps);
        try
            globalF = thisF((this_f-15):(this_f+30));
        catch
            globalF = [];
        end
        all_trials = [all_trials;globalF];
    end
    
%     all_trials = all_trials(1:55,:);
    
    spine_fig = figure, plot(all_trials'), hold on, plot(mean(all_trials),'Color','k','LineWidth',4)
    yy = ylim;
    line([16 16],[yy])
    line([23 23],[yy])
    title(['spine ' num2str(handpicked(sp))])
    set(gcf,'Position',[200 200 600 1200])
    saveas(spine_fig,[saveDir_allDendrites2 filesep 'D' num2str(this_dendrite) '_spine' num2str(handpicked(sp)) '_posthocSTA.fig']);
    
end
% 
% for sp = 1:size(spinesF,2)
%     
%     thisF = spinesF(:,sp)';
%     
%     figure, plot(thisF), hold on
%     scatter(ps_times, repelem(2,numel(ps_times)))
%     title(num2str(handpicked(sp)))
% 
% end

%% new section: post hoc check if visual responses
sync_saveDir2 = ['Z:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) '\sigTargets'];
sync_posthoc = extractPaq_dendrites_targetPermPOSTHOC_vis(animal, session, this_dendrite, 1, nFrames, sync_saveDir2, saveDir2)

allSP=readtable(['E:\Data\M' num2str(animal) '\S' num2str(session) '\D' num2str(this_dendrite) 'globalF.csv']);
thisF=table2array(allSP(:,[2]))';
figure, plot(thisF)%hold on; plot(thisF2); figure,plot(thisF-thisF2);thisF1=thisF;thisF=thisF-thisF2;

ps_times = [];
vis_times = [];
vis_order = [];
for block = 1:numel(sync_posthoc.sync_blocks)
    these = sync_posthoc.sync_blocks{block}.framesWithPhotostim;
    ps_times = [ps_times (these+(block-1)*8000)]
    
    these = sync_posthoc.sync_blocks{block}.vis_frames;
    vis_times = [vis_times (these+(block-1)*8000)]
    
    these = sync_posthoc.sync_blocks{block}.angle_by_trial;
    vis_order = [vis_order these]

end


% downsample to 4x avg?
% ps_times = floor(ps_times/4);
% vis_times = floor(vis_times/4);

% ps_times = floor(ps_times/4);
% vis_times = floor(vis_times/4);


figure, 
angles = unique(vis_order);
for angle = 1:numel(unique(vis_order))
    this_angle = angles(angle);
    trials = find(vis_order == this_angle);
    times = vis_times(trials);
    
    STAs_this_angle = [];
    for t = times
        this_STA = thisF([t-30:t+60]);
        STAs_this_angle = [STAs_this_angle;this_STA];
    end
    subplot(1,numel(unique(vis_order)), angle)
    plot((STAs_this_angle')), hold on
    plot(mean(STAs_this_angle),'Color','k','LineWidth',3)
    yy = ylim;
    line([31 31],[yy],'Color','b')
    line([60 60],[yy],'Color','b')
    
    line([36 36],[yy],'Color','r')
    line([43 43],[yy],'Color','r')
    title(num2str(this_angle))
end

STAs = [];
for ps_trial = ps_times
    this_STA = thisF([ps_trial-30:ps_trial+60]);
    STAs = [STAs;this_STA];
end
figure, 
plot((STAs')), hold on
plot(mean(STAs),'Color','k','LineWidth',3)
line([31 31],[yy],'Color','b')
line([38 38],[yy],'Color','b')



%% clear everything specific to this dendrite: 
clear this_dendrite this_soma_location thisROI WF_soma_px soma_AE pairwise_distances




%% new section DH - 2021-05-17: sparse noise sync extraction 
outputDir=saveDir2;

this_dendrite=2
sync_saveDir =['Z:\Data\M1098\S4\D' num2str(this_dendrite) '\sigTargets'];
outputDir=['E:\Data\M1098\S4\targetPerm_sigTargets\D' num2str(this_dendrite)];
if ~exist(outputDir)
    mkdir(outputDir)
end
getSparseNoiseSync(animal, session, this_dendrite, sync_saveDir, outputDir)

% soma SN 

animal = 1110;
session = 3;
this_soma = 1;

sync_saveDir =['Z:\Data\M' num2str(animal) '\S' num2str(session)];
outputDir=['E:\Data\M' num2str(animal) '\S' num2str(session) '\SN_soma'];

if ~exist(outputDir)
    mkdir(outputDir)
end

% regular sparse noise 
getSparseNoiseSync_vSoma(animal, session, this_soma, sync_saveDir, outputDir)

% sparse noise with different noise files for each block 

getSparseNoiseSync_vSomaMultipleNoiseFiles(animal, session, this_soma, sync_saveDir, outputDir)

% 4x4 grid of small gratings to localize RF 
expt =3;
nplanes = 1;
getGratingSync_grid(animal, session, expt, nplanes,this_soma) % should mostly work but need to adapt to correct format of x and y shifts 

% Full Field grating before small grating grid 
expt = 1;
nplanes = 2;
getGratingSync_FFtuning(animal, session, expt, nplanes,this_soma) % should mostly work but need to adapt to correct format of x and y shifts 





