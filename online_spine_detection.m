    
% animal=1076;
% session=3;
dendrite=this_dendrite;

%SELECT ALL IN ROI MANAGER >> MORE >> MULTI-MEASURE. CHECK BOTH BOXES. SAVE
cd(['E:\Data/M' num2str(animal) '/S' num2str(session) ])
allSP=readtable(['D' num2str(dendrite) 'allSP2.csv']);
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

allSP=readtable(['D' num2str(dendrite) 'allSP3.csv']);
this=table2array(allSP(:,[3:4:size(allSP,2)]));

spinesF2=NaN(size(this,1),numSpines); 
for spine=1:numSpines
    thisTrace=this(:,spine);
    thisF0=running_percentile(thisTrace,1350,10);%thisF0=thisF0';
    thisdF=thisTrace-thisF0; thisdFF=thisdF./thisF0;
    spinesF2(:,spine) = thisdFF;
    fprintf('spine %0d F obtained\n',spine)
end
spinesF=[spinesF; spinesF2];

figure,plot(1:size(spinesF,1),spinesF)

figure, hold on 
for sp = 1:size(spinesF,2)
    thisF = smooth(spinesF(:,sp));
    plot(thisF+((sp-1)*2));
end
lims = ylim(gca);
ylim([-1 lims(2)])

%% GET SYNC MAKE RASTER

thisSync=sync_this_dendrite;

blockSize=8000;numBlocks=3;
stimTimes=[];blankTimes=[];
for block=1:numBlocks
    stimFrames = thisSync.sync_blocks{block}.frames_ps;
    stimFrames = stimFrames + (blockSize*(block-1));
    stimTimes=[stimTimes stimFrames];
    blankFrames = thisSync.sync_blocks{block}.frames_blank;
    blankFrames = blankFrames + (blockSize*(block-1));
    blankTimes=[blankTimes blankFrames];
end

allStimOrder=thisSync.stim_parameters_config.final_order(:,1:numBlocks);
stimOrder = reshape(allStimOrder,1,[]);
rasterDuration = 4;halfDur=rasterDuration/2;
approxFR=round(15);%approxFR=round(15/4);
stimRasterF = NaN(numSpines,length(stimTimes),(approxFR*rasterDuration)+1);
blankRasterF = NaN(numSpines,length(blankTimes),(approxFR*rasterDuration)+1);

for spine=1:numSpines
    thisdFF=spinesF(:,spine)';
    for st=1:length(blankTimes)
        t=blankTimes(st);
        if t<(approxFR*(halfDur))
            starts = 1;
            ends = t+(approxFR*halfDur);pad = NaN(1,((approxFR*rasterDuration)+1)-ends);
            blankRasterF(spine,st,:)=[pad thisdFF(starts:ends)];
        else 
            starts = t-(approxFR*halfDur);
            ends = t+(approxFR*halfDur);
            if ends>length(thisdFF)
                pad=NaN(1,ends-length(thisdFF));
                blankRasterF(spine,st,:)=[thisdFF(starts:end) pad];
            else
                blankRasterF(spine,st,:)=thisdFF(starts:ends);
            end
        end
    end
    for st=1:length(stimTimes)
        t=stimTimes(st);       
        if t<(approxFR*(halfDur))
            starts = 1;
            ends = t+(approxFR*halfDur);pad = NaN(1,((approxFR*rasterDuration)+1)-ends);
            stimRasterF(spine,st,:)=[pad thisdFF(starts:ends)];
        else 
            starts = t-(approxFR*halfDur);
            ends = t+(approxFR*halfDur);
            if ends>length(thisdFF)
                pad=NaN(1,ends-length(thisdFF));
                stimRasterF(spine,st,:)=[thisdFF(starts:end) pad];
            else
                stimRasterF(spine,st,:)=thisdFF(starts:ends);
            end
        end
    end
end
meanRast=squeeze(mean(stimRasterF(:,:,:),1));

figure,subplot(2,1,1)
plot(1:size(meanRast,2),meanRast,'-');hold on
line([9 9],[0 3],'color','r','LineWidth',2)
line([11 11],[0 3],'color','r','LineWidth',2)
plot(mean(meanRast,1),'color','k','LineWidth',2)
title(['all stims M' num2str(animal) 'S' num2str(session) 'D' num2str(dendrite)])

subplot(2,1,2)
plot(mean(meanRast,1),'color','k','LineWidth',2)
% savefig(['allStimsFig'])
% savefig(['/Users/mehmetfisek/Desktop/darkness_nonClust_registered/' ...
%     ' M' num2str(animal) 'S' num2str(session) 'D' num2str(dendrite) 'allStimsFig'])

%% GET TARGETS

% allTargets=horzcat(thisSync.stim_parameters_config.target_groups(:).these);
% allTargets=allTargets';

allTargets=vertcat(thisSync.stim_parameters_config.target_groups(:).these);
allCoords=vertcat(thisSync.stim_parameters_config.target_groups(:).these_coordinatesXYZ);
targCoord=[];
for t=1:length(unique(allTargets))
    targCoord=[targCoord; allCoords(find(allTargets==t,1),:)];
end

% allTargetGroups=vertcat(thisSync.stim_parameters_config.target_groups(:).these);
allTargetGroups=[];
for j=1:length(thisSync.stim_parameters_config.target_groups(:))
    this=thisSync.stim_parameters_config.target_groups(j).these;
    if length(this)<n_targets
        this2=[this; NaN(n_targets-length(this),1)];
    else
        this2=this;
    end
    allTargetGroups=[allTargetGroups this2];
end
   allTargetGroups=allTargetGroups';
    
[sorted,sortind]=sort(stimOrder);

% trying pseudo-targets
figure,scatter(targCoord(:,1),-1*targCoord(:,2));hold on
scatter(targCoord(:,1)+2,-1*targCoord(:,2)-1.5)
% pseudoTargets=[targCoord(:,1)+2 targCoord(:,2)-2 targCoord(:,3)];
pseudoTargets=[targCoord(:,1)+2 targCoord(:,2)-1.5 targCoord(:,3)];

numPlanes=length(unique(targCoord(:,3)));
plane_spacing;
plane_length=zeros(1,numPlanes);
pseudoPairs=[];
for tplane=1:numPlanes
    thisPlaneROIs=find(targCoord(:,3)==plane_spacing(tplane));pseudoPair=[];
    plane_length(tplane)=length(thisPlaneROIs);
    dists=squareform(pdist([targCoord(thisPlaneROIs,:);pseudoTargets(thisPlaneROIs,:)]));
    dists=dists(1:length(thisPlaneROIs),length(thisPlaneROIs)+1:end);dists=(dists<3);
    for p=1:length(thisPlaneROIs)
            these=find(dists(p,:));
            pseudoPair(p,1:length(these))=these+sum(plane_length(1:tplane-1));
    end
    pseudoPairs=[pseudoPairs; pseudoPair];
end
    
    
%% get mean spine signal for global events and detect hit targers
% respStart=21;respEnd=25;
% baseStart=16;baseEnd=20;
respStart=39;respEnd=45; %was to 13
% respStart=36;respEnd=45; %was to 13
% respStart=32;respEnd=38; %was to 13
baseStart=22;baseEnd=30;
respsMn=mean(meanRast(:,respStart:respEnd),2);
basesMn=mean(meanRast(:,baseStart:baseEnd),2);
deltsMn=respsMn-basesMn;


%% GET EFFECTIVE GLOBAL TARGETS

bythr=[];
for thr=1:4

   globalSDthresholds=[1 1.5 2 10];% multiple of SD to set as threshold above mean for picking trials
%      globalSDthresholds=[3 3.5 4 4.5]
%     globalReliabilities=[0.1 0.2];

    globalSDthreshold=globalSDthresholds(thr);% multiple of SD to set as threshold above mean for picking trials
    globalReliability=0.10;

    trialGrpMn=find(deltsMn>(mean(deltsMn)+(globalSDthreshold*std(deltsMn)))); %
    stimGrpMn=stimOrder(trialGrpMn);
%     targetsMn=vertcat(thisSync.stim_parameters_config.target_groups(stimGrpMn).these);
    targetsMn=allTargetGroups(stimGrpMn,:);   
    
    targetsMncount=histcounts(targetsMn,'BinEdges',0.5:1:length(targCoord)+0.5);
    numtrialsMn=size(targetsMn,1)
    ratioThresMn=max(2,round(numtrialsMn*globalReliability));
    hitTargetsMn=find(targetsMncount>ratioThresMn);
    hitCoordsMn=targCoord(hitTargetsMn,:);
    bythr=[bythr; [thr numel(hitTargetsMn)]];
end
bythr


    thr=4;
    globalSDthreshold=globalSDthresholds(thr);% multiple of SD to set as threshold above mean for picking trials
    globalReliability=0.5;

    trialGrpMn=find(deltsMn>(mean(deltsMn)+(globalSDthreshold*std(deltsMn)))); %
    stimGrpMn=stimOrder(trialGrpMn);
%     targetsMn=vertcat(thisSync.stim_parameters_config.target_groups(stimGrpMn).these);
    targetsMn=allTargetGroups(stimGrpMn,:);   
    targetsMncount=histcounts(targetsMn,'BinEdges',0.5:1:length(targCoord)+0.5);
    numtrialsMn=size(targetsMn,1);ratioThresMn=max(2,round(numtrialsMn*globalReliability));
    hitTargetsMn=find(targetsMncount>ratioThresMn);
    hitCoordsMn=targCoord(hitTargetsMn,:);

%% GET EFFECTIVE SPINE TARGETS

spineReliability=0.10; %proportion of trials passing threshold that must contain hit   
spineSDthresholds=[0.5 1 1.5 2];
spineSDthresholds=[2 2.5 3 3.5];

bythr=[];hits1={};
for thr=1:4
    spineSDthreshold=spineSDthresholds(thr);   
    for sp=1:numSpines
        thisRast=squeeze(stimRasterF(sp,:,:));
        resps=mean(thisRast(:,respStart:respEnd),2);
        bases=mean(thisRast(:,baseStart:baseEnd),2);
        baseSDs=std(thisRast(:,baseStart:baseEnd),[],2);
        delts=resps-bases;
        deltsZ=resps-bases;
        mnThres=mean(deltsMn)+(1*std(deltsMn));
        loneThres=mean(delts(find(deltsMn<mnThres)))+...
            (spineSDthreshold*std(delts(find(deltsMn<mnThres))));
        trialGrp1=find(delts>loneThres & deltsMn<mnThres);
        stimGrp1=stimOrder(trialGrp1);
        
%         targets1=vertcat(thisSync.stim_parameters_config.target_groups(stimGrp1).these);
        targets1=allTargetGroups(stimGrp1,:); 
        
        targets1count=histcounts(targets1,'BinEdges',0.5:1:length(targCoord)+0.5);
        numtrials1=size(targets1,1);%numtrials2=size(targets2,1);
        ratioThres1=max(2,round(numtrials1*spineReliability));%ratioThres2=max(2,round(numtrials2*0.1));
        hitTargets1=find(targets1count>ratioThres1);hitCoords1=targCoord(hitTargets1,:);
        hits1{sp}=hitTargets1;
        [r1,c1]=ind2sub(size(allTargetGroups),find(ismember(allTargetGroups(:),hits1{sp})));
        r1=sortind(r1);
        [r,c]=ind2sub(size(allTargetGroups),find(ismember(allTargetGroups(:),hitTargetsMn)));
        rm=sortind(r);
    end
    bythr=[bythr; [thr numel(unique([hits1{:}])) numel(unique([hits1{:}]))/numSpines]];
end

bythr

    hits1={};
    thr=2
    spineSDthreshold=spineSDthresholds(thr);
%     spineReliability=0.2; %proportion of trials passing threshold that must contain hit   

    for sp=1:numSpines
        thisRast=squeeze(stimRasterF(sp,:,:));
        resps=mean(thisRast(:,respStart:respEnd),2);
        bases=mean(thisRast(:,baseStart:baseEnd),2);
        delts=resps-bases;
        mnThres=mean(deltsMn)+(1*std(deltsMn));
        loneThres=mean(delts(find(deltsMn<mnThres)))+...
            (spineSDthreshold*std(delts(find(deltsMn<mnThres))));

        trialGrp1=find(delts>loneThres & deltsMn<mnThres);
        stimGrp1=stimOrder(trialGrp1);
        
%         targets1=vertcat(thisSync.stim_parameters_config.target_groups(stimGrp1).these);
        targets1=allTargetGroups(stimGrp1,:); 
        
        targets1count=histcounts(targets1,'BinEdges',0.5:1:length(targCoord)+0.5);
        numtrials1=size(targets1,1);%numtrials2=size(targets2,1);
        ratioThres1=max(2,round(numtrials1*spineReliability));%ratioThres2=max(2,round(numtrials2*0.1));
        hitTargets1=find(targets1count>ratioThres1);
        
%         hitTargets1=hitTargets1(3);
        
        hitCoords1=targCoord(hitTargets1,:);
        hits1{sp}=hitTargets1;
        [r1,c1]=ind2sub(size(allTargetGroups),find(ismember(allTargetGroups(:),hits1{sp})));  
        r1=sortind(r1);        
        r1hotTrials=r1(ismember(r1,trialGrp1));
        [r,c]=ind2sub(size(allTargetGroups),find(ismember(allTargetGroups(:),hitTargetsMn)));
        rm=sortind(r);
        
        if ~isempty(r1)
            figure, 
            subplot(1,3,1)
            plot(1:61,squeeze(stimRasterF(sp,1:20:size(stimRasterF,2),:)),'color',[0.7 0.7 0.7],'LineWidth',1);hold on

            plot(1:61,squeeze(stimRasterF(sp,r1,:)),'color',[250 110 50]/255,'LineWidth',2);
            plot(1:61,squeeze(stimRasterF(sp,r1hotTrials,:)),'color',[170 0 40]/255,'LineWidth',2);
            plot(1:61,mean(squeeze(stimRasterF(sp,r1,:)),1),'color','b','LineWidth',2);
        
            line([30 30],[0 3],'color','g','LineWidth',2)
            line([38 38],[0 3],'color','g','LineWidth',2)
            xlim([0 60])
            subplot(1,3,2)
            plot(deltsMn,delts,'ob','MarkerSize',5);hold on
            plot(deltsMn(r1),delts(r1),'o','color',[250 110 50]/255,'MarkerSize',9)
            plot(deltsMn(r1hotTrials),delts(r1hotTrials),'o','color',[255 0 40]/255,'MarkerSize',9)
            plot(deltsMn(rm),delts(rm),'*g','MarkerSize',7)
    %         xlim([-1 4]);ylim([-2 4])
            line([mnThres mnThres],[-2 4],'color',[0.5 0.5 0.5])
            line([-1 4],[loneThres loneThres],'color',[0.5 0.5 0.5])
            line([-1 4],[-1 5],'color',[0.5 0.5 0.5])
            set(gcf,'Position',[100 300 2400 800])   
            subplot(1,3,3)
            view(2);hold on
            scatter3(hitCoords1(:,1),hitCoords1(:,2),hitCoords1(:,3),48,'or');
            xlim([0 512]);ylim([0 512])
            scatter3(hitCoordsMn(:,1),hitCoordsMn(:,2),hitCoordsMn(:,3),48,'*g');
            title(['spine ' num2str(sp) ' active min ' num2str(ratioThres1) ' trials out of ' num2str(numtrials1)])
        end
    end

%% COLLECT


handpicked = [12];

storing = [spineReliability*ones(length(handpicked),1) ...
    spineSDthresholds(thr)*ones(length(handpicked),1) ...
    handpicked' ];
hpHitnum=[];
for h=1:length(handpicked)
    hpHitnum=[hpHitnum; length(hits1{handpicked(h)})];
end
storing=[storing hpHitnum];
% stored=[];
stored=[stored;storing];

% handPcoords=[];
for s=1:length(handpicked)
    h=handpicked(s);n=length(hits1{h});
    this=[h*ones(n,1) hits1{h}'  targCoord(hits1{h},:)];
    handPcoords=[handPcoords;this];
end



globalHits=[];%hitTargetsMn;
globalHitCoords=[];%hitCoordsMn;

% spineHits=unique([hits1{:}]);
% 
% spineHits=[spineHits hits1{handpicked}];
% spineHits=[hits1{handpicked}];

spineHits=unique(handPcoords(:,2)');


spineHitCoords=targCoord(spineHits,:);
[numel(globalHits) numel(spineHits)]
newHits=[globalHits spineHits];
newCoords=[globalHitCoords; spineHitCoords]
length(spineHits)/numSpines
numTarget=length(newHits)
groupSize = numTarget;

% groupSize=min(12,round(numTarget*0.5))

numGroups=1; % 

% run this for CONFIRMATION blocks and ANN if you want to stim just a
% single location 
enrichedTargets = [];
for g=1:numGroups
    theseTs=randperm(numTarget,groupSize);
    enrichedTargets(g).group=g;
    enrichedTargets(g).these=theseTs;
    enrichedTargets(g).these_coordinatesXYZ=newCoords(theseTs,:);
end
% 
% nums=[nums; numel(globalHits) numel(spineHits) numel(spineHits)/numSpines]

save(['D' num2str(dendrite) 'enrichedTargets'],'enrichedTargets','numTarget')
    
save(['D' num2str(dendrite) 'enrichedTargets_all'], 'enrichedTargets','globalHits','globalHitCoords',...
    'spineHits','spineHitCoords','handpicked','stored','handPcoords');
    
save(['D' num2str(dendrite) 'handPcoords'],'handPcoords')
    

% run this to stim a COLUMN around sig targets 

spineHits=unique(handPcoords(1,2)');

spineHits=[845]



spineHitCoords=targCoord(spineHits,:);
[numel(globalHits) numel(spineHits)]
newHits=[globalHits spineHits];
newCoords=[globalHitCoords; spineHitCoords]
length(spineHits)/numSpines
numTarget=length(newHits)
groupSize = numTarget;

% groupSize=min(12,round(numTarget*0.5))

numGroups=1; % 

% run this for CONFIRMATION blocks and ANN if you want to stim just a
% single location 
enrichedTargets = [];
for g=1:numGroups
    theseTs=randperm(numTarget,groupSize);
    enrichedTargets(g).group=g;
    enrichedTargets(g).these=theseTs;
    enrichedTargets(g).these_coordinatesXYZ=newCoords(theseTs,:);
end


px2um = 3.75;
dist_fromCenterTarget = 9; % in um
axial_offset = 30;
x_offset = round(cos(deg2rad(30))*dist_fromCenterTarget/px2um);
y1_offset = round(dist_fromCenterTarget/px2um);
y2_offset = round(sin(deg2rad(30))*dist_fromCenterTarget/px2um);

do3D = 0;

enrichedTargets = [];
for g=1:numGroups
    theseTs=randperm(numTarget,groupSize);
    enrichedTargets(g).group=g;
    enrichedTargets(g).these=theseTs;
    
    mainTarget_px = newCoords(theseTs,:);
    
    these_coordinatesXYZ = [];

    for t = 1:size(mainTarget_px,1)
        these_ccords = mainTarget_px(t,:);
        % compute xy coords for prism stimulation pattern in same plane
        pt_in0{t}= [these_ccords(1) these_ccords(2) these_ccords(3)];
        pt_in1{t}= [these_ccords(1) these_ccords(2)-y1_offset these_ccords(3)];
        pt_in2{t}= [these_ccords(1)-x_offset these_ccords(2)+y2_offset these_ccords(3)];
        pt_in3{t}= [these_ccords(1)+x_offset these_ccords(2)+y2_offset these_ccords(3)];
        
        these_coordinatesXYZ = [these_coordinatesXYZ;[pt_in0{t};pt_in1{t};pt_in2{t};pt_in3{t}]];

        if do3D
            % compute xy coords for prism stimulation pattern in above and
            % below planes
            pt_ab0{t}= [these_ccords(1) these_ccords(2) these_ccords(3)+axial_offset];
            pt_ab1{t}= [these_ccords(1) these_ccords(2)+y1_offset these_ccords(3)+axial_offset];
            pt_ab2{t}= [these_ccords(1)-x_offset these_ccords(2)-y2_offset these_ccords(3)+axial_offset];
            pt_ab3{t}= [these_ccords(1)+x_offset these_ccords(2)-y2_offset these_ccords(3)+axial_offset];
            
            pt_bl0{t}= [these_ccords(1) these_ccords(2) these_ccords(3)-axial_offset];
            pt_bl1{t}= [these_ccords(1) these_ccords(2)+y1_offset these_ccords(3)-axial_offset];
            pt_bl2{t}= [these_ccords(1)-x_offset these_ccords(2)-y2_offset these_ccords(3)-axial_offset];
            pt_bl3{t}= [these_ccords(1)+x_offset these_ccords(2)-y2_offset these_ccords(3)-axial_offset];
            

            
            these_coordinatesXYZ = [these_coordinatesXYZ;...
                [pt_ab0{t};pt_ab1{t};pt_ab2{t};pt_ab3{t};...
                pt_bl0{t};pt_bl1{t};pt_bl2{t};pt_bl3{t}]];

        end

    end
    enrichedTargets(g).these_coordinatesXYZ=these_coordinatesXYZ;
end
figure, subplot(121),
scatter3(these_coordinatesXYZ(1,1),these_coordinatesXYZ(1,2),these_coordinatesXYZ(1,3),50,'filled','r'), hold on
scatter3(these_coordinatesXYZ(2:end,1),these_coordinatesXYZ(2:end,2),these_coordinatesXYZ(2:end,3),50,'filled','b')

spiralSize = 16;
spiralSize_px = spiralSize/px2um;
subplot(122)
for trg = 1:size(these_coordinatesXYZ,1)
    p = nsidedpoly(1000, 'Center', [these_coordinatesXYZ(trg,1) these_coordinatesXYZ(trg,2)], 'Radius', spiralSize_px/2);
    plot(p, 'FaceColor', 'r'), hold on
end

save(['D' num2str(dendrite) 'enrichedTargets_column'],'enrichedTargets','numTarget')
    
save(['D' num2str(dendrite) 'enrichedTargets_all_column'], 'enrichedTargets','globalHits','globalHitCoords',...
    'spineHits','spineHitCoords','handpicked','stored','handPcoords');
    