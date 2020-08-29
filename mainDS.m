%% Main DS script
% created 8/28/2018 by Alex Tiriac
% last update: 8/28/2020
%
% Version history:
% v3: Sends some data to another folder that will be used by a table
% v2: Fixed two bugs
%   -When finding pref close to 360, it used to snag to 315, this is fixed
%   -The polarPlots were plotting maxdF instead of vector sum, also fixed
%   -added ability to remove "bad cells"
% v1.1: plotting pref direction with vector sum now
% v1: does not discriminate b/w ON or OFF
%
%
%%% Table of contents:
%%%
%%% Inputs
%%% User loading of files
%%% Give an index for ON, ON-OFF, and BAD cells
%%% Corrects directions based on pre calibration (see evernote exp file)
%%% use ROI mask to find z profile at each cell
%%% Find max intensity at every bar presentation (not discriminating between ON and OFF)
%%% Tuning curves
%%% Plot DSI v VS   and     maxDF histogram to check threshold values
%%% map out where the DS cells are in the retina
%%% Plot DS polar plots
%%% Plot raw traces
%%% Save ROIint and textFile as a variable (and sends it to a folder)
%
%


%% Inputs

shouldIGetRidofBackF = 0; % Set to 0 if you don't want to background subtract light stim
frameRate = 1.48;  % 1.48  
DSI_thresh = 0;
vecSum_thresh = 0;
maxDF_thresh = 0;

lengthOfFOV = 850/2; %For 2x zoom, it's 850um/2

OnOffCells =  [2,14,15,24,29,30,31,34,35,41,43,46];
OnCells = [4,5,6,32,36,42,44,45];
badCells = [1,3,7:13,16:23,25,26,27,28,33,37,38,39,40];

goodCells = [OnOffCells,OnCells];

VS_axis = 0.5;





%%Direction correction from DMD comp to rig
%
%To do this correction:
%1) Look up photo of MOM dmd calib
%2) draw retina with actual dir on paper oriented in the rig (ventral=270)
%3) look at 0 on postit note and assign correct value of oriented retina.
%note: real0 refers to the direction of stim on retina in rig
%more note, left side is postit, right side is retina
%




real0 = 135;
real45 = 90;
real90 = 45;
real135 = 0;
real180 = 315;
real225 = 270;
real270 = 225;
real315 = 180;
%End dir correction




%% User loading of files
%Loading the calcium movie
disp('pick the calcium imaging file'); 
[movie_name, path_name] = uigetfile([path, 'experiment', '\*.tif']);

%Loading the mask
disp('pick the calcium ROI mask file'); 
[roi_mask_file, path_name] = uigetfile([path, 'experiment', '\*.tif']);

%Loading the text files
disp('pick the text file with DS directions')
[textF, path_name_text] = uigetfile([path, 'experiment', '\*.txt']);

textFile = readtable(textF);
textFileArray = table2array(textFile);
textFileArray = [textFileArray, textFileArray(:,8)-textFileArray(:,7)];
StimDuration = mean(textFileArray(:,9));
[num_trials dummy] = size(textFileArray);

%Loading the movie
n_frames = numel(imfinfo(movie_name));
[height, width] = size(imread(movie_name,1));
movieMat = zeros(height, width, n_frames);
for k = 1:n_frames
    movieMat(:,:,k) = imread(movie_name, k);
end
%End loading the movie

%Loading the mask
roi_mask = imread(roi_mask_file);
countmask=bwlabel(roi_mask);
%End loading the mask

disp('Done loading movie'); %Status update

%% Give an index for ON, ON-OFF, and BAD cells

num_cells = max(max(countmask));

cellID = ones(num_cells,1);

cellID(OnOffCells) = 2;
cellID(OnCells) = 3;
cellID(badCells) = 4;


%% Correct directions
newInd = zeros(24,1);
replaceDir = [real0; real45; real90; real135; real180; real225; real270; real315];
temp = unique(textFileArray(:,1));

for i = 1:numel(temp)
    indices = textFileArray(:,1) == temp(i);
    newInd(indices) = replaceDir(i);
end
textFileArrayOld = textFileArray;
textFileArray(:,1) = newInd; 


%% use ROI mask to find z profile at each cell

for i=1:n_frames
    partprops=regionprops(countmask,movieMat(:,:,i),'Area','MeanIntensity','Centroid'); 
    for j = 1:max(max(countmask))
        roiInt(j,i) = partprops(j).MeanIntensity;
        cellCenters(j,:) = partprops(j).Centroid; %not centered and in pixels
%         cellCenters(j,:) = partprops(j).Centroid/256*lengthOfFOV-0.5*lengthOfFOV; %centered and in um
    end
    
end



[num_cells, num_frames] = size(roiInt);





%% Find max intensity at every bar presentation (not discriminating between ON and OFF)

timeOfStim = mean(textFileArray(:,8)-textFileArray(:,7));
timeOfStim_frames = ceil(timeOfStim * frameRate);
winWvf = 3;

barResp = zeros(num_cells,num_trials);
wvf_resp = zeros(num_cells, timeOfStim_frames+2*winWvf ,num_trials); % waveforms of all responses
wvf_resp_allComb = [];

for i = 1:num_trials
    
    stimFrameStart = floor(textFileArray(i,7)*frameRate);
    stimFrameEnd = ceil(textFileArray(i,8)*frameRate);
    barResp(:,i) = max(roiInt(:,stimFrameStart:stimFrameEnd),[],2);
    
    wvf_resp(:,:,i) = roiInt(:,stimFrameStart-winWvf:stimFrameStart+timeOfStim_frames+winWvf-1);
    wvf_resp_allComb = [wvf_resp_allComb; roiInt(:,stimFrameStart-winWvf:stimFrameStart+timeOfStim_frames+winWvf-1)];
    
end

ind0 = find(textFileArray(:,1) == 0);
ind45 = find(textFileArray(:,1) == 45);
ind90 = find(textFileArray(:,1) == 90);
ind135 = find(textFileArray(:,1) == 135);
ind180 = find(textFileArray(:,1) == 180);
ind225 = find(textFileArray(:,1) == 225);
ind270 = find(textFileArray(:,1) == 270);
ind315 = find(textFileArray(:,1) == 315);



%% Tuning curves
disp('Computing DS'); %Status update

%Zeros start
DSI = zeros(num_cells,1);
vecSum = zeros(num_cells,1);
vecTheta = zeros(num_cells,1);
maxDF = zeros(num_cells,1);

rhos_all = zeros(num_cells, 9);
rhos1_all = zeros(num_cells, 9);
rhos2_all = zeros(num_cells, 9);
rhos3_all = zeros(num_cells, 9);
rhos_var = zeros(num_cells, 9);
%Zeroes end




thetas = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 0]; %Theta values for 0 45 90 135 180 225 270 315


for i = 1:num_cells
    
    rhos = [mean(barResp(i,ind0)), mean(barResp(i,ind45)), mean(barResp(i,ind90)), mean(barResp(i,ind135)), mean(barResp(i,ind180)), mean(barResp(i,ind225)), mean(barResp(i,ind270)), mean(barResp(i,ind315)), mean(barResp(i,ind0))];
    rhos1 = [barResp(i,ind0(1)), barResp(i,ind45(1)), barResp(i,ind90(1)), barResp(i,ind135(1)), barResp(i,ind180(1)), barResp(i,ind225(1)), barResp(i,ind270(1)), barResp(i,ind315(1)), barResp(i,ind0(1))];
    rhos2 = [barResp(i,ind0(2)), barResp(i,ind45(2)), barResp(i,ind90(2)), barResp(i,ind135(2)), barResp(i,ind180(2)), barResp(i,ind225(2)), barResp(i,ind270(2)), barResp(i,ind315(2)), barResp(i,ind0(2))];
    rhos3 = [barResp(i,ind0(3)), barResp(i,ind45(3)), barResp(i,ind90(3)), barResp(i,ind135(3)), barResp(i,ind180(3)), barResp(i,ind225(3)), barResp(i,ind270(3)), barResp(i,ind315(3)), barResp(i,ind0(3))];
    rhos_var(i,:) = [var(barResp(i,ind0)), var(barResp(i,ind45)), var(barResp(i,ind90)), var(barResp(i,ind135)), var(barResp(i,ind180)), var(barResp(i,ind225)), var(barResp(i,ind270)), var(barResp(i,ind315)), var(barResp(i,ind0))];
  
    rhosNorm = rhos/sum(rhos(1:8));
    [x,y] = pol2cart(thetas,rhosNorm);
    vectSumX = sum(x(1:8));
    vectSumY = sum(y(1:8));
    [the, rho] = cart2pol(vectSumX,vectSumY);
    if the <0
        the = 2*pi+the;
    end
    
    %what is the pref dir?
    
    if the > 5.8905
        prefInd = 1;
    else
        [temp,prefInd] = min(abs(thetas-the));
    end
    prefDir = thetas(prefInd);
    nullInd = prefInd - 4;
    if nullInd < 1
        nullInd = 8+nullInd;
    end
    
    DSI_temp = (rhos(prefInd)-rhos(nullInd))/(rhos(prefInd)+rhos(nullInd));
    
    
    %Store for cells
    DSI(i,1) = DSI_temp;
    vecSum(i,1) = rho;
    vecTheta(i,1) = the;
    maxDF(i,1) = max(rhos); %What was the max DF
    
    rhos_all(i,:) = rhos;
    rhos1_all(i,:) = rhos1;
    rhos2_all(i,:) = rhos2;
    rhos3_all(i,:) = rhos3;  
    %End store for cells
    

end

%% Plot DSI v VS   and     maxDF histogram to check threshold values

disp('generating plots'); %Status update


figure('Name', 'cell DS, VS, & max dF', 'position', [100, 100, 1000, 400])
subplot(1,2,1)
plot(vecSum(goodCells,1), DSI(goodCells,1),'k.');
axis([0 1 0 1])
xlabel('Norm Vector Sum')
ylabel('DSI')
subplot(1,2,2)
histogram(maxDF(goodCells), 'BinLimits', [0, 3], 'BinWidth', 0.1)




%% Map out where the DS cells are in the retina
% This is where I use the DSI and VS and DF thresholds specified at start

ind = DSI>DSI_thresh & vecSum>vecSum_thresh & maxDF > maxDF_thresh;

% Map out Vector angle regardless of cell type
colorDS = countmask;
for i = 1:length(ind)
    colorDS(colorDS == i) = vecTheta(i);
end
% end of map out vector angle regardless of cell type


% Now make masks of countmask for ON-OFF and ON cells
cellID_DS = cellID .* ind;

indOnOffDS = find(cellID_DS == 2);
indOnDS = find(cellID_DS == 3);
indAllDS = [indOnOffDS; indOnDS];

whereDS = countmask;
whereDS_OnOff = countmask;
whereDS_On = countmask;

for i = 1:length(indAllDS)
    whereDS(whereDS == indAllDS(i)) = 0.5;
end

whereDS(whereDS ~= 0.5) = 0;
whereDS(whereDS == 0.5) = 1;

for i = 1:length(indOnOffDS)
    whereDS_OnOff(whereDS_OnOff == indOnOffDS(i)) = 0.5;
end

whereDS_OnOff(whereDS_OnOff ~= 0.5) = 0;
whereDS_OnOff(whereDS_OnOff == 0.5) = 1;

for i = 1:length(indOnDS)
    whereDS_On(whereDS_On == indOnDS(i)) = 0.5;
end

whereDS_On(whereDS_On ~= 0.5) = 0;
whereDS_On(whereDS_On == 0.5) = 1;

% End of make masks of countmask for ON-OFF and ON cells

whereDS = whereDS.*colorDS;
whereDS_OnOff = whereDS_OnOff.*colorDS; %Where the OnOff cells are w/ pref angle
whereDS_On = whereDS_On.*colorDS; %Where the On cells are w/ pref angle


%Now plot it all

figure('Name', 'DS cell location', 'position', [0, 44, 1200, 300])
s1 = subplot(1,3,1);
imagesc(whereDS) 
cmap = hsv;
cmap(1,:) = 0;
colormap(cmap);
caxis([0 2*pi]);
s2 = subplot(1,3,2);
imagesc(whereDS_OnOff) 
cmap = hsv;
cmap(1,:) = 0;
colormap(cmap);
caxis([0 2*pi]);
s3 = subplot(1,3,3);
imagesc(whereDS_On) 
cmap = hsv;
cmap(1,:) = 0;
colormap(cmap);
caxis([0 2*pi]);

%% Plot DS polar plots

figure('Name', 'DS polar plots', 'position', [0, 44, 1200, 300])
s1 = subplot(1,3,1);
polarplot([vecTheta(indAllDS) vecTheta(indAllDS)]',[zeros(length(indAllDS),1) vecSum(indAllDS)]','r');
rlim([0 VS_axis]);
set(gca,'thetaticklabel',{[]})
set(gca,'rticklabel',{[]})
thetaticks([])
rticks([])
s2 = subplot(1,3,2);
polarplot([vecTheta(indOnOffDS) vecTheta(indOnOffDS)]',[zeros(length(indOnOffDS),1) vecSum(indOnOffDS)]','r');
rlim([0 VS_axis]);
set(gca,'thetaticklabel',{[]})
set(gca,'rticklabel',{[]})
thetaticks([])
rticks([])
s3 = subplot(1,3,3);
polarplot([vecTheta(indOnDS) vecTheta(indOnDS)]',[zeros(length(indOnDS),1) vecSum(indOnDS)]','r');
rlim([0 VS_axis]);
set(gca,'thetaticklabel',{[]})
set(gca,'rticklabel',{[]})
thetaticks([])
rticks([])

%% Plot DS polar plots

nBins = 24;
rLimVar = 15;

figure('Name', 'DS polar histograms', 'position', [0, 44, 1200, 300])
s1 = subplot(1,3,1);
polarhistogram(vecTheta(indAllDS),nBins);
set(gca,'thetaticklabel',{[]})
set(gca,'rticklabel',{[]})
thetaticks([])
rticks([])
% rlim([0 rLimVar])
s2 = subplot(1,3,2);
polarhistogram(vecTheta(indOnOffDS),nBins);
set(gca,'thetaticklabel',{[]})
set(gca,'rticklabel',{[]})
thetaticks([])
rticks([])
% rlim([0 rLimVar])
s3 = subplot(1,3,3);
polarhistogram(vecTheta(indOnDS),nBins);
set(gca,'thetaticklabel',{[]})
set(gca,'rticklabel',{[]})
thetaticks([])
rticks([])
% rlim([0 rLimVar])


%% Plot raw traces

reordered_ind = [ind0 ; ind45 ; ind90 ; ind135 ; ind180 ; ind225 ; ind270 ; ind315];
wvf_resp_reordered = wvf_resp(:,:,reordered_ind);
wvf_resp_mean = [mean(wvf_resp_reordered(:,:,1:3),3),mean(wvf_resp_reordered(:,:,4:6),3),mean(wvf_resp_reordered(:,:,7:9),3),mean(wvf_resp_reordered(:,:,10:12),3),mean(wvf_resp_reordered(:,:,13:15),3),mean(wvf_resp_reordered(:,:,16:18),3),mean(wvf_resp_reordered(:,:,19:21),3),mean(wvf_resp_reordered(:,:,22:24),3)];

xaxisMax = 160;

%For ON-OFFs
[~,idx] = sort(vecTheta(indOnOffDS));
sortedOnOffs = [vecTheta(indOnOffDS), indOnOffDS];
sortedOnOffs = sortedOnOffs(idx,:);

figure('Name', 'OnOff traces', 'Units', 'Normalized', 'position', [0, 0.1, 0.25, 0.025*length(indOnOffDS)])
hold
for i = 1:length(indOnOffDS)
    plot(wvf_resp_mean(sortedOnOffs(i,2),:)+i, 'k', 'LineWidth',1)
end
ylim([0, length(indOnOffDS)+1])
xlim([0, xaxisMax]);
xticks([10:20:150])
xticklabels({'0' '45' '90' '135' '180' '225' '270' '315'})

%For ONs
[~,idx2] = sort(vecTheta(indOnDS));
sortedOns = [vecTheta(indOnDS), indOnDS];
sortedOns = sortedOns(idx2,:);

figure('Name', 'On traces', 'Units', 'Normalized', 'position', [0.25, 0.1, 0.25,0.025*length(indOnDS)])
hold
for i = 1:length(indOnDS)
    plot(wvf_resp_mean(sortedOns(i,2),:)+i, 'k', 'LineWidth',1)
end
ylim([0, length(indOnDS)+1])
xlim([0, xaxisMax]);
xticks([10:20:150])
xticklabels({'0' '45' '90' '135' '180' '225' '270' '315'})


figure, imagesc(countmask)

[x,y,idxOrdered, C] = classDS(vecTheta,vecSum, 4);


%% Save ROIint and textFile as a variable


% pathName = 'C:\Users\Olive\Dropbox\DS project\3_compiledResults\Data\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'cellCenters', 'textFileArray')

% pathName = 'C:\Users\Olive\Dropbox\DS project\8_Satb2\SummaryData\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'textFileArray')


% pathName = 'C:\Users\Olive\Dropbox\DS project\1_Data\190612 piezo2 gfp\PiezoResultsFolder\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'textFileArray')

% pathName = 'C:\Users\Olive\Dropbox\DS project\Data 2P Ca2+ imaging for Alex\P13-14\b140815 (retina2)\retina # 2\analysis\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'textFileArray')


% pathName = 'C:\Users\Olive\Dropbox\DS project\6_singularityMapping\190816 vglut2creG6s mapping singularities\dataPool\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'cellCenters', 'textFileArray')

% pathName = 'C:\Users\Olive\Dropbox\DS project\13_Hb9Cx36KO\summaryData\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'cellCenters', 'textFileArray')

% pathName = 'C:\Users\Olive\Dropbox\DS project\3_compiledResults\Data\b2ko\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'cellCenters', 'textFileArray')

% pathName = 'C:\Users\Olive\Dropbox\DS project\3_compiledResults\Data\DHBE exp\';
% MovieNameNoExt = movie_name(1:end-4);
% save([pathName,MovieNameNoExt,'.mat'], 'roiInt', 'cellCenters', 'textFileArray')




%% old useless code

% whereDS = countmask;
% vecTheta_colors = vecTheta(indAllDS);
% vecSum_DS = vecSum(indAllDS);
% 
% for i = 1:length(ind)  %First make a mask of DS cells "whereDS"
%     whereDS(whereDS == ind(i)) = 0.5;
% end
% 
% whereDS(whereDS ~= 0.5) = 0;
% whereDS(whereDS == 0.5) = 1;
% 
% 
% colorDS = countmask.*whereDS;  %DS cells now colored based on pref dir
% for i = 1:length(ind)
%     colorDS(indAllDS) = vecTheta_colors(i);
% end
% 
% figure('Name', 'DS cell location', 'position', [0, 44, 650, 600])
% s1 = subplot(2,2,1);
% imagesc(countmask)  %All cells
% hold
% s2 = subplot(2,2,2);
% imagesc(colorDS)    %DS cells color-coded based on pref dir
% cmap = hsv;
% cmap(1,:) = 0;
% colormap(cmap);
% caxis([0 2*pi]);
% colormap(s1,parula);
% s3 = subplot(2,2,3);
% polarplot([vecTheta_colors vecTheta_colors]',[zeros(length(vecTheta_colors),1) vecSum_DS]','r');
% rlim([0 VS_axis]);
% s4 = subplot(2,2,4);
% polarplot([vecTheta_colors vecTheta_colors]',[zeros(length(vecTheta_colors),1) ones(length(vecTheta_colors),1)]','r');
% 
% figure('Name', 'ClusterPlot', 'position', [0, 44, 400, 400])
% polarplot([vecTheta_colors vecTheta_colors]',[zeros(length(vecTheta_colors),1) vecSum_DS]','r');
% rlim([0 VS_axis]);
% set(gca,'thetaticklabel',{[]})
% set(gca,'rticklabel',{[]})
% thetaticks([])
% rticks([])



%% functions!

function [x,y,idxOrdered, C] = classDS(prefDir,vecSum, numClusters)

    [x,y] = pol2cart(prefDir,vecSum);
    [xNorm,yNorm] = pol2cart(prefDir,ones(length(prefDir),1));

    [idx,C] = kmeans([xNorm yNorm],numClusters);

    idxOrdered = nan(length(idx),1);

    % Reorder grps based on max X or Y
    %Temporal = 1
    [dummy,id] = max(C(:,1));
    idxOrdered(idx == id) = 1;
    %Dorsal = 2
    [dummy,id] = max(C(:,2));
    idxOrdered(idx == id) = 2;
    %Nasal = 3
    [dummy,id] = min(C(:,1));
    idxOrdered(idx == id) = 3;
    %Ventral = 4
    [dummy,id] = min(C(:,2));
    idxOrdered(idx == id) = 4;
    %Done with Reorder

end









