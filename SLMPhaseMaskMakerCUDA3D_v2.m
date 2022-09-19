function varargout = SLMPhaseMaskMakerCUDA3D_v2(varargin)
% Lloyd Russell 2015
% Implements a weighted GS algorithm using Martin Persson's HOTlab DLL
% (https://github.com/MartinPersson/HOTlab)

% Parameters
% =========================================================================
% None         : open a load file GUI
% TargetsImage : 512*512 array or { arrays }
% FileName     : string or { strings }
% Points       : n*4 array (x,y,z,I) or { arrays }
% all_Galvo_Positions: n*2 phase mask group galvo positions 
% -------------------------------------------------------------------------
% Save         : boolean (default=true)
% SaveName     : string or { strings }
%                    N.B. must include 'tif' extension  --- should fix this
% SaveFormat   : string, 'tiff' or 'bmp' (default='tiff')
% BitDepth     : integer, 8 or 16 (default=8)
% -------------------------------------------------------------------------
% DoTransform  : boolean (default=true)
% TransformFile    : string
% -------------------------------------------------------------------------
% Iterations   : integer
% -------------------------------------------------------------------------
% SliceSpacing : integer, step size between slices of input 3D stack
% FocalSlice   : integer, which slice of 3D stack is the normal focal plane
% ManualZ      : integer or 1d array of length n spots or { array }

% Returns
% =========================================================================
% PhaseMasks
% TransformedSLMTargets
% -------------------------------------------------------------------------
% Saved files:
% |- PhaseMasks                     < results folder
% |----- (filename)_CUDAphase.ext   < the computed phase mask
% |----- InputTargets               < folder of input targets
% |----- TransformedTargets         < folder of transformed targets

% NB: -um puts spots deeper in sample than native focal plan of objective

p = inputParser;
p.addParameter('TargetsImage', []);
p.addParameter('SaveName', '');
p.addParameter('SaveDirectory', '');
p.addParameter('FileName', '');
p.addParameter('Do2DTransform', 0); % was true 
p.addParameter('Do3DTransform', 1); % was true 
p.addParameter('TransformFile', '');
p.addParameter('TransformFile3D', '');
p.addParameter('TransformFile3D_9GPs', ''); % DH
p.addParameter('Save', true);
p.addParameter('SaveFormat', 'tiff');
p.addParameter('BitDepth', '8');
p.addParameter('Iterations', 100); % was 100, had to change it for making large phase masks, can be set higher for up to 1024x1024 
p.addParameter('Is3D', false);
p.addParameter('SliceSpacing', 1);
p.addParameter('FocalSlice', 0);
p.addParameter('ManualZ', []);
p.addParameter('Points',[]);
p.addParameter('AutoAdjustWeights',1);
p.addParameter('SteepnessFudgeFactor',[]);
p.addParameter('all_Galvo_Positions', []),... % DH

parse(p, varargin{:});

FocalSlice{1} = p.Results.FocalSlice;
SliceSpacing = p.Results.SliceSpacing;

% read yaml file to get settings
yaml = ReadYaml('settings.yml');

% Load HOTlab DLL
HologramLibraryName = yaml.HologramLibraryName;
if ~libisloaded(HologramLibraryName)
    loadlibrary([HologramLibraryName '.dll'])
end

% decide which transform path to use (2D)
Do2DTransform = p.Results.Do2DTransform;
if isempty(p.Results.TransformFile)
    TransformFile = yaml.TransformFile;
else
    TransformFile = p.Results.TransformFile;
end
if Do2DTransform
    load(TransformFile, 'tform');  % loads 'tform'
end


% improved calibration with 9 galvo positions

Do3DTransform = p.Results.Do3DTransform;
if isempty(p.Results.TransformFile3D)
    TransformFile3D_9GP = yaml.TransformFile3D_9GPs;
else
    TransformFile3D_9GP = p.Results.TransformFile3D_9GPs;
end
if Do3DTransform
    load(TransformFile3D_9GP);  % loads 'T'
end

% decide which transform path to use (3D)
% Do3DTransform = p.Results.Do3DTransform;
% if isempty(p.Results.TransformFile3D)
%     TransformFile3D = yaml.TransformFile3D;
% else
%     TransformFile3D = p.Results.TransformFile3D;
% end
% if Do3DTransform
%     load(TransformFile3D);  % loads 'T'
% end


% see if AutoAdjustWeights passed as input, if not take value from yaml file
if isempty(p.Results.AutoAdjustWeights)
    AutoAdjustWeights = yaml.AutoAdjustWeights;
else
    AutoAdjustWeights = p.Results.AutoAdjustWeights;
end
% AutoAdjustWeights = 1
if AutoAdjustWeights
    load(['C:\Users\User\Dropbox\Bruker3\SLM\Weighting' ...
        filesep yaml.WeightingFile], 'W');  % loads 'W'
    disp(['autoadjusting weights based on the following file: ' ['C:\Users\User\Dropbox\Bruker3\SLM\Weighting' ...
        filesep yaml.WeightingFile]])
end

% fudge weights steepness (slope of fit)
if isempty(p.Results.SteepnessFudgeFactor)
    SteepnessFudgeFactor = yaml.SteepnessFudgeFactor;
else
    SteepnessFudgeFactor = p.Results.SteepnessFudgeFactor;
end

% see if BitDepth passed as input, if not take value from yaml file
if isempty(p.Results.BitDepth)
    BitDepth = yaml.SLM_BitDepth;
else
    BitDepth = p.Results.BitDepth;
end



% INPUT: Manual points
% -------------------------------------------------------------------------
if ~any(strcmpi(p.UsingDefaults, 'Points'))
    if any(strcmpi(p.UsingDefaults, 'SaveName')) && p.Results.Save
        warning('Must provide save name when passing in raw points')
        return
    end
    
    if iscell(p.Results.Points)
        Points = p.Results.Points;
        NumFiles = length(p.Results.Points);
        SaveName = p.Results.SaveName;
    else
        NumFiles = 1;
        Points{1} = p.Results.Points;
        SaveName{1} = p.Results.SaveName;
    end
    
    for idx = 1:NumFiles
        num_targets = numel(Points{idx}(:,1));
        num_planes = ceil(range(unique(Points{idx}(:,3)))+1);
        InputTargets{idx} = zeros(512,512,num_planes);
        for i = 1:num_targets
            row = round(Points{idx}(i,2));
            col = round(Points{idx}(i,1));
            slice = round(Points{idx}(i,3) - min(Points{idx}(:,3)) + 1);
            val = Points{idx}(i,4);
            InputTargets{idx}(row,col,slice) = val;
        end
        FocalSlice{idx} = -min(Points{idx}(:,3)) + 1;
    end
    
    
% INPUT: Image array provided
% -------------------------------------------------------------------------
elseif ~any(strcmpi(p.UsingDefaults, 'TargetsImage'))
    if ~any(strcmpi(p.UsingDefaults, 'SaveName'))  % save name must also be provided
        if iscell(p.Results.TargetsImage)  % passed multiple images
            InputTargets = p.Results.TargetsImage;
            SaveName = p.Results.SaveName;
        else
            NumFiles = 1;
            InputTargets{1} = p.Results.TargetsImage;
            SaveName{1} = p.Results.SaveName;
        end
    else
        warning('Must provide save name when passing in image data')
    end
    
    
% INPUT: Filepath provided
% -------------------------------------------------------------------------
elseif ~any(strcmpi(p.UsingDefaults, 'FileName'))
    FileName = p.Results.FileName;
    if iscell(FileName)  % passed multiple filenames
        NumFiles = length(FileName);
        InputTargets = cell(NumFiles,1);
        
        for idx = 1:length(FileName)
            info = imfinfo(FileName{idx});
            num_images = numel(info);
            InputTargets{idx} = zeros(512,512,num_images);
            for k = 1:num_images
                InputTargets{idx}(:,:,k) = imread(FileName{idx}, k);
            end
        end
    else
        NumFiles = 1;
        info = imfinfo(FileName);
        num_images = numel(info);
        InputTargets{1} = zeros(512,512,num_images);
        for k = 1:num_images
            InputTargets{1}(:,:,k) = imread(FileName, k);
        end
        FileName = {FileName};
    end
    SaveName = FileName;
    
    
% INPUT: none, open load file dialog
% -------------------------------------------------------------------------
elseif any(strcmpi(p.UsingDefaults, 'TargetsImage')) && any(strcmpi(p.UsingDefaults, 'FileName'))
    [FileName, PathName] = uigetfile('*.tif*', 'Select the targets file(s)', 'MultiSelect', 'on');
    if PathName == 0
        warning('Cancelled')
        return
    end
    transform_choice = questdlg('Transform points?', 'User input required', 'Yes', 'No', 'Yes');
    if strcmpi(transform_choice, 'No')
        Do2DTransform = false;
    end
    
    % load images into array
    cd(PathName);
    InputTargets = {};
    if iscell(FileName)  % selected multiple files
        NumFiles = length(FileName);
        InputTargets = cell(NumFiles,1);
        for idx = 1:length(FileName)
            info = imfinfo(FileName{idx});
            num_images = numel(info);
            InputTargets{idx} = zeros(512,512,num_images);
            for k = 1:num_images
                InputTargets{idx}(:,:,k) = imread(FileName{idx}, k);
                FocalSlice{idx} = 1;  %fix this later!
            end
        end
    else
        NumFiles = 1;
        info = imfinfo(FileName);
        num_images = numel(info);
        InputTargets{1} = zeros(512,512,num_images);
        for k = 1:num_images
            InputTargets{1}(:,:,k) = imread(FileName, k);
        end
        FileName = {FileName};
    end
    SaveName = FileName;
end

% rebuild input targets if provided manual Z coordinates
if ~any(strcmpi(p.UsingDefaults, 'ManualZ'))
    if iscell(p.Results.ManualZ)
        manualZ = p.Results.ManualZ;
    else
        manualZ{1} =  p.Results.ManualZ;
    end
    if numel(manualZ) < NumFiles
        manualZ = repmat(manualZ,1,NumFiles);
    end
    
    for idx = 1:NumFiles
        % Find target coordinates
        [y, x, z] = ind2sub(size(InputTargets{idx}),find(InputTargets{idx}));
        
        % get target intensities
        [~,~,I] = find(InputTargets{idx});
        
        % replace Z with manual Z
        if numel(manualZ{idx}) ~= numel(z)
            manualZ{idx} = repmat(manualZ{idx},1,numel(z));
        end
        z = manualZ{idx};
        
        % build new input image
        num_targets = numel(x);
        num_planes = range(unique(z))+1;
        InputTargets{idx} = zeros(512,512,num_planes);
        for i = 1:num_targets
            row = y(i);
            col = x(i);
            slice = z(i) - min(z) + 1;
            val = I(i);
            InputTargets{idx}(row,col,slice) = val;
        end
        FocalSlice{idx} = -min(z) + 1;
    end
end


SLMsize = [2048 2048];  % x y


% Start CUDA
deviceId    = 0;
h_pSLMstart = zeros(SLMsize);  % the starting phase mask
LUT         = [];
% cuda_error1 = calllib(HologramLibraryName,'startCUDA', h_pSLMstart, deviceId);


% Process all points/files (do transforms, make phase masks and save images)
PhaseMasks = cell(NumFiles,1);
TransformedSLMTargets = cell(NumFiles,1);
for f = 1:NumFiles
%     tic;
    disp(['computing phase mask for group ' num2str(f) ' of ' num2str(NumFiles)])

    % Find target coordinates
    [y_i, x_i, z_i] = ind2sub(size(InputTargets{f}),find(InputTargets{f}));
    
    % get target intensities
    [~,~,I] = find(InputTargets{f});
    
    % adjust z coordinates
    z_i = (z_i - FocalSlice{f}) * SliceSpacing;
    
    % get how many different planes
    im_dims = size(InputTargets{f});
    num_dims = length(im_dims);
    if num_dims > 2
        num_planes_input = im_dims(3);
    else
        num_planes_input = 1;
    end
    
    % transform coordinates from imaging space to SLM space
    if Do3DTransform
        points = [x_i y_i z_i];
    % note DH: have to give it all Ts for all Gavlo position to be able to
    % interpolate 
        
%         if ~ exist(p.Results.all_Galvo_Positions(1))
%             warning('no galvo position provided, setting to 256, 256 by default')
%             transformedPoints = applyAffineTransform3D_GalvoOffset(points, [256,256], T); % 
% 
%         else
           transformedPoints = applyAffineTransform3D_GalvoOffset(points, p.Results.all_Galvo_Positions(f,:), T); % 
%         end 
        
%         keyboard
        x = transformedPoints(:,1);
        y = transformedPoints(:,2);
        z = transformedPoints(:,3);
%         keyboard
    else
        if Do2DTransform
            [x, y] = transformPointsForward(tform, x_i, y_i);
        else
            x = x_i; y = y_i;
        end
        z = z_i;
    end
    
    % Weight the target pixel intensity by distance from zero order
    if AutoAdjustWeights % modified DH 2019-08-21
        
       distances = pairwiseDistance([x-256 y-256], [0 0]); % is this correct? may have to be center of SLM zero order which is 700 550 or sth 
%         distances = pairwiseDistance([points(:,1) points(:,1)], [256 256]); %DH 2019-08-21


        
        % fudge weights steepness (slope of fit)
        W.p_edited = W.p;
%         W.p_edited(1) = W.p_edited(1) * SteepnessFudgeFactor; % original Lloyd's code 
        W.p_edited(1) = W.p_edited(1); % edited DH to implement fudge factor for 2nd order fit 
%         SteepnessFudgeFactor
        % estimate intensites of spots based on calibration data
        estimatedIntensity = polyval(W.p_edited, distances);
        MAX = polyval(W.p, 0);
        normEstimatedIntensity = estimatedIntensity / MAX; % percent of power at hypothetical zero order position 
        
        % compute weights from estimated intensites
%         weights = (1 ./ normEstimatedIntensity);  % original line from LLoyd, DH implemented linear subtraction for new SLM... subtract to ensure linear correction. division does not.
        % DH compute weights 2020-01-07
        offsets = -1*(normEstimatedIntensity-max(normEstimatedIntensity))*SteepnessFudgeFactor;
        weights = 1+offsets;
        
        % apply weights
        I_weighted = I .* weights;
        
        if any(I_weighted<=0)
            disp('ERROR: Some spot weights <= 0')
            I_weighted(I_weighted<=0) = min(I_weighted(I_weighted>0));
        end
        
        % save weights
        I_orig = I;
        I = I_weighted;
%         figure, scatter3(x,y,I)
    end
    
    % Check for and remove targets outside of SLM space
%     if any(x<1 | x>512 | y<1 | y>512)
%         warning('At least one transformed target spot is outside of addressable SLM space. The offending spots will be removed.')
%         idx_remove = (x<1 | x>512 | y<1 | y>512);
%         x(idx_remove) = [];
%         y(idx_remove) = [];
%         z(idx_remove) = [];
%         I(idx_remove) = [];
%     end
    
    
    
    % MAKE PHASE MASK
    SLMsize         = [yaml.SLM_Pixels_X yaml.SLM_Pixels_Y];  % x y
    pSLMstart       = zeros(SLMsize);  % the starting phase mask
    pSLM            = zeros(SLMsize);
    x_spots         = x - 256;  %(SLMsize(1)/2);  % subtract 256 because centre of image is 0,0 (not 256,256)
    y_spots         = y - 256;  %(SLMsize(2)/2);  % subtract 256 because centre of image is 0,0 (not 256,256)
    z_spots         = z;
    I_spots         = I;  % intensity comes from input image
    N_spots         = length(x_spots);
    calcIntensities = 0;  % unused
    calcTime        = 0;  % always 0
    method          = 1;
    N_iterations    = 8; % was 20 
    useGpu          = 1;
    max_n_spots     = 100;

    errCode1 = calllib(HologramLibraryName, 'Create_generator_cl',...
        SLMsize(1),SLMsize(2),max_n_spots,N_iterations,useGpu);

    [errCode2, ~, ~, ~, ~, ~, ~, pSLM, ~, ~] = calllib(HologramLibraryName,'Generate_hologram_cl',...
        N_spots, x_spots, y_spots, z_spots, I_spots, N_iterations, method, pSLMstart, pSLM, calcIntensities, calcTime);

    errCode3 = calllib(HologramLibraryName, 'Destroy_generator_cl');



%      figure; imagesc(pSLM); axis equal; axis tight; box off

 
    
    % Convert
    pSLM = uint8(pSLM);  
    PhaseMasks{f} = pSLM(449:1600, 65:1984); % crop to rectangular size 
    
%     figure; imagesc(PhaseMasks{f}); axis equal; axis tight; box off

%     
%     % Generate Hologram parameters
%     h_test         = [];
%     h_test_START   = [];
%     h_pSLM         = zeros(SLMsize);
%     x_spots        = x - 256; %(SLMsize(1)/2);  % subtract 256 because centre of image is 0,0 (not 256,256)
%     y_spots        = y - 256; %(SLMsize(2)/2);  % subtract 256 because centre of image is 0,0 (not 256,256)
%     z_spots        = z;
%     I_spots        = I;  % intensity comes from input image
%     N_spots        = length(x_spots);
%     N_iterations   = p.Results.Iterations;
%     h_obtainedAmps = [];
%     h_obtainedAmps_START = [];
% 
%     if any(z_spots ~= 0)
%         method = 1;
%     else
%         method = 2;
%     end
%     disp(num2str(x_spots))
%     disp(num2str(y_spots))
%     disp(num2str(z_spots))
%     
%     method = 0;
%     N_iterations = 1;
%     
%     
% %     N_iterations
% %     method
%     % 0: Complex addition of "Lenses and Prisms", no optimization (3D)
%     % 1: Weighted Gerchberg-Saxton algorithm using Fresnel propagation (3D)
%     % 2: Weighted Gerchberg-Saxton algorithm using fast fourier transforms (2D)
%     
%     % Apply corrections
%     % [cuda_error, h_AberrationCorr, h_LUTPolCoeff, h_LUT_uc] = calllib(LibraryName, 'Corrections',...
%     %     UseAberrationCorr, h_AberrationCorr, UseLUTPol, PolOrder, h_LUTPolCoeff, saveAmplitudes, alpha, DCborderWidth, UseLUT, h_LUT_uc);
%     
%     % Generate Hologram
%     
%     
%     
%     
%     cuda_error1 = calllib(HologramLibraryName, 'startCUDA', h_pSLMstart, deviceId); % added DH to re-initialize before every phase mask
%     
%     [cuda_error2,~,h_pSLM,~,~,~,~,h_obtainedAmps] = calllib(HologramLibraryName,'GenerateHologram',...
%         h_test_START, h_pSLMstart, x_spots, y_spots, z_spots, I_spots, N_spots, N_iterations, h_obtainedAmps_START, method);
%          figure; imagesc(h_pSLM)
% 
%       
% %      calllib(HologramLibraryName, 'Generate_hologram_cl', unsigned int n_spots, x_positions, y_positions, z_positions, ...
% %          intensities, n_iterations, method, starting_phases, hologram_image, calc_intensities, ...
% %          calc_time_us )
% %     
% 
% 
%     n_spots = N_spots;
%     x_positions = x_spots;
%     y_positions = y_spots;
%     z_positions = z_spots;
%     intensities = I_spots;
%     n_iterations = N_iterations;
%     method = method;
%     starting_phases = h_pSLMstart;
%     hologram_image = [];
%     calc_intensities = [];
%     calc_time_us = [];
% 
%      calllib(HologramLibraryName, 'Generate_hologram_cl', n_spots, x_positions, y_positions, z_positions, ...
%          intensities, n_iterations, method, starting_phases, hologram_image, calc_intensities, ...
%          calc_time_us )
%     
%     
%     % crop for rectangular SLM
% %     h_pSLM = h_pSLM(1:1152, 1:1920);
%     
%     h_pSLM = h_pSLM(449:1600, 65:1984);
% 
%     
    
    
    
%     % Convert
%     if BitDepth == 8
%         PhaseMasks{f} = uint8(h_pSLM / 256);  % convert to 8 bit
%     else
%         PhaseMasks{f} = uint16(h_pSLM);
%     end
    
   %     figure; imagesc(h_pSLM)
%     title(['iterations = ' num2str(N_iterations) ', method: ' num2str(method)])
    
    % Build transformed targets image
    x = round(x);
    y = round(y);
    z = round(z);  % neccessary to round z units?
    
    % get how many different planes
    num_planes_output = range(unique(z))+1;
    
    
%         TransformedSLMTargets{f} = zeros(512, 512, num_planes_output);
%         for i = 1:length(x)
%             y(y < 1) = 1; % DH: set coordinates outside image to 0 
%             x(x < 1) = 1;
%             TransformedSLMTargets{f}(y(i), x(i), z(i)-min(z)+1) = I(i);  
%         end
% 
%     
    
    % Save
    if p.Results.Save
        SaveDirectory = p.Results.SaveDirectory;
        if ~isempty(SaveDirectory)
            StartingDirectory = pwd;
            cd(SaveDirectory)
        end
        
        % Make output directories
        if ~exist([pwd filesep 'PhaseMasks'], 'dir')
            mkdir([SaveDirectory filesep 'PhaseMasks']);
        end
        if ~exist([pwd filesep 'PhaseMasks' filesep 'InputTargets'], 'dir')
            mkdir([pwd filesep 'PhaseMasks' filesep 'InputTargets']);
        end
%         if ~exist([pwd filesep 'PhaseMasks' filesep 'TransformedTargets'], 'dir')
%             mkdir([pwd filesep 'PhaseMasks' filesep 'TransformedTargets']);
%         end
        
        if isempty(strfind(SaveName{f},'.tiff')) && ~isempty(strfind(SaveName{f},'.tif'))
            imwrite(PhaseMasks{f}, ['PhaseMasks' filesep strrep(SaveName{f}, '.tif', ['_CUDAphase.' p.Results.SaveFormat])]);
        elseif ~isempty(strfind(SaveName{f},'.tiff'))
            imwrite(PhaseMasks{f}, ['PhaseMasks' filesep strrep(SaveName{f}, '.tiff', ['_CUDAphase.' p.Results.SaveFormat])]);
        else
            sprintf('--- Warning: unknown image format')
        end
%         for k = 1:num_planes_input
%             if k == 1
%                 writemode = 'overwrite';
%             else
%                 writemode = 'append';
%             end
%             imwrite(InputTargets{f}(:,:,k), ['PhaseMasks' filesep 'InputTargets' filesep strrep(SaveName{f}, '.tif', '_InputTargets.tif')], 'writemode', writemode);
%         end
        
%         for k = 1:num_planes_output
%             if k == 1
%                 writemode = 'overwrite';
%             else
%                 writemode = 'append';
%             end
%             imwrite(TransformedSLMTargets{f}(:,:,k), ['PhaseMasks' filesep 'TransformedTargets' filesep strrep(SaveName{f}, '.tif', '_TransformedTargets.tif')], 'writemode', writemode);
%         end
%         
        % move back to starting directory
        if ~isempty(SaveDirectory)
            cd(StartingDirectory)
        end
    end
%     time_stop = toc;
    %     disp(['[' num2str(f) ' of ' num2str(NumFiles) '] Done in ' num2str(round(time_stop*1000*1000)/1000) ' ms'])
end

% % Stop CUDA
% cuda_error3 = calllib(HologramLibraryName,'stopCUDA');
% unloadlibrary(HologramLibraryName);  % tidy up, unload dll

% Outputs if requested
if nargout
    varargout{1} = PhaseMasks;
    varargout{2} = TransformedSLMTargets;
end
