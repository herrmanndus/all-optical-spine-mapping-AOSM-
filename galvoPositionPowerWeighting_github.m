function galvo_scale_factor = galvoPositionPowerWeighting_github(g_x, g_y)
% Dustin Herrmann 2020-06-19
% weight power by the galvo offset
% extreme gavlo positions lead to power dropoff, compensate for that here 


%% calibration and create lookup table file 
% 
% galvoGridMargin = 50; % was 100 
% galvoGridNumPoints = 5; % was 3
% [galvoX, galvoY] = meshgrid(linspace(galvoGridMargin, 512-galvoGridMargin,galvoGridNumPoints),linspace(galvoGridMargin, 512-galvoGridMargin,galvoGridNumPoints)) ;  % in pixels
% numGalvoLocations = numel(galvoX);
% 
% isSpiral = 'True';
% spiralDiameter = 85;
% MarkPoints_GPLMaker(round(galvoX(:)), round(galvoY(:)), isSpiral, spiralDiameter, 3, ['Transform3DGalvoPositions_' timestamp]);
% 
% 
% % measure powers at all 25 positions. Note, the order is first rows then
% % columns 
% 
% galvoX = galvoX(:);
% galvoY = galvoY(:);
% 
% measurements = ...
%     [150 183 192 183 155 ...
%     200 215 222 217 196 ...
%     212 220 230 221 212 ...
%     208 220 223 220 207 ...
%     178 205 204 206 180];
% 
% cftool
% % do simple interpolant fit 

%% apply 

try
    load('C:\Users\User\Dropbox\Bruker3\PowerUtilities\powerWeightingGalvos.mat')
catch
    load('C:\Users\Dustin\Dropbox\Bruker3\PowerUtilities\powerWeightingGalvos.mat')
end

maxi = fittedmodel_galvoP(256,256);
current = fittedmodel_galvoP(g_x,g_y);
galvo_scale_factor = maxi/current;


end 