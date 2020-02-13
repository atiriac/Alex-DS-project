
%% Inputs

frameRate = 1.48;  % 1.48  
DSI_thresh = 0.3;
vecSum_thresh = 0;
maxDF_thresh = 0;

lengthOfFOV = 850/2; %For 2x zoom, it's 850um/2

OnOffCells =  [1,2,3,9,10,11,12,14,16,18,19,23,27,30,32,34,35,36,37,38,40,41,42,43,44,53,55,56,57,58,59,61,62,65,66,67,68,69,70,71,73,74,81,82,85,88,94,95,98,105,106,107,109,111,112,114,115,117,118,128,129,138,139,141,142,145];
OnCells = [4,13,17,20,24,46,51,52,60,83,86,87,89,90,102,116,119,120,121,123,124,125,126,127,130,131,132,134,140];
badCells = [5,6,7,8,15,21,22,25,26,28,29,31,33,39,45,47,48,49,50,54,63,64,72,75,76,77,78,79,80,84,91,92,93,96,97,99,100,101,103,104,108,110,113,133,135,136,137,143,144];

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


% real0 = 90;
% real45 = 135;
% real90 = 180;
% real135 = 225;
% real180 = 270;
% real225 = 315;
% real270 = 0;
% real315 = 45;

real0 = 0;
real45 = 45;
real90 = 90;
real135 = 135;
real180 = 180;
real225 = 225;
real270 = 270;
real315 = 315;
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

cellCenters(:,2) = cellCenters(:,2) * (-1);


% get rid of bad cells
roiInt(badCells,:) = [];
cellCenters(badCells,:) = [];

[num_cells, num_frames] = size(roiInt);


%% calc DS
[DSI,vecSum, vecTheta, wvf_resp_mean, sumVar, rhos_all] = calculateDS(roiInt, textFileArray, frameRate);


% determine DSI of cells for randomized dir using bootstrapping
shuffledDSI = nan(length(DSI),1000);
parfor j = 1:1000
    shuffledText = textFileArray;
    shuffledText(:,1) = shuffleTrialDirections(textFileArray(:,1),8, 3);

    [shuffledDSI(:,j), dummy, dummy, dummy, dummy, dummy] = calculateDS(roiInt, shuffledText, frameRate);
end

%Determine DSI significance by comparing real DSI to randomize DSI
for i = 1:num_cells
    
    [val,idx]=min(abs(DSI(i)-sort(shuffledDSI(i,:))));
    DSIsig(i,1) = idx/1000;
end

DSIsigThresh = 0.95;

% indDS = find(DSIsig > 0.95);

% Get rid of non DS cells
% DSI(DSIsig < DSIsigThresh) = [];
% vecTheta(DSIsig < DSIsigThresh) = [];
% vecSum(DSIsig < DSIsigThresh) = [];


%% Functional classification of DS cells based on preferred directions

silValue = silTest(vecTheta);
figure, plot(mean(silValue))

numClusters = 4;
[x,y,idxClass, C] = funClassDS(vecTheta,vecSum, numClusters);

figure('Name', 'Functional Clusters', 'Position', [300 400 1000 400])
subplot(1,2,1)
scatter(x,y,[],idxClass)
hold on

subplot(1,2,2)
polarplot([vecTheta(idxClass == 1) vecTheta(idxClass == 1)]',[zeros(length(vecTheta(idxClass == 1)),1) vecSum(idxClass == 1)]','m');
hold on
polarplot([vecTheta(idxClass == 2) vecTheta(idxClass == 2)]',[zeros(length(vecTheta(idxClass == 2)),1) vecSum(idxClass == 2)]','b');
polarplot([vecTheta(idxClass == 3) vecTheta(idxClass == 3)]',[zeros(length(vecTheta(idxClass == 3)),1) vecSum(idxClass == 3)]','g');
polarplot([vecTheta(idxClass == 4) vecTheta(idxClass == 4)]',[zeros(length(vecTheta(idxClass == 4)),1) vecSum(idxClass == 4)]','y');

%% vector flow fields

xQuiver = cellCenters(:,1);
yQuiver = cellCenters(:,2);
[uQuiver vQuiver] = pol2cart(vecTheta,vecSum);

classNum = 1;

figure('Name', 'quiver', 'Position', [400 400 500 420])
quiver(xQuiver(idxClass == classNum),yQuiver(idxClass == classNum),uQuiver(idxClass == classNum),vQuiver(idxClass == classNum));
axis([0 255 -255 0])

X = xQuiver(idxClass == classNum);
Y = yQuiver(idxClass == classNum);
[Xq,Yq] = meshgrid(1:1:250, -250:1:-1);
Uq = griddata(X,Y,uQuiver(idxClass == classNum),Xq,Yq);
Vq = griddata(X,Y,vQuiver(idxClass == classNum),Xq,Yq);

figure('Name', 'flowField', 'Position', [400 400 500 420])
streamslice(Xq, Yq, Uq, Vq)
axis([0 255 -255 0])


% 
% %% Find max intensity at every bar presentation (not discriminating between ON and OFF)
% 
% timeOfStim = mean(textFileArray(:,8)-textFileArray(:,7));
% timeOfStim_frames = ceil(timeOfStim * frameRate);
% winWvf = 3;
% 
% barResp = zeros(num_cells,num_trials);
% wvf_resp = zeros(num_cells, timeOfStim_frames+2*winWvf ,num_trials); % waveforms of all responses
% wvf_resp_allComb = [];
% 
% for i = 1:num_trials
%     
%     stimFrameStart = floor(textFileArray(i,7)*frameRate);
%     stimFrameEnd = ceil(textFileArray(i,8)*frameRate);
%     barResp(:,i) = max(roiInt(:,stimFrameStart:stimFrameEnd),[],2);
%     
%     wvf_resp(:,:,i) = roiInt(:,stimFrameStart-winWvf:stimFrameStart+timeOfStim_frames+winWvf-1);
%     wvf_resp_allComb = [wvf_resp_allComb; roiInt(:,stimFrameStart-winWvf:stimFrameStart+timeOfStim_frames+winWvf-1)];
%     
% end
% 
% ind0 = find(textFileArray(:,1) == 0);
% ind45 = find(textFileArray(:,1) == 45);
% ind90 = find(textFileArray(:,1) == 90);
% ind135 = find(textFileArray(:,1) == 135);
% ind180 = find(textFileArray(:,1) == 180);
% ind225 = find(textFileArray(:,1) == 225);
% ind270 = find(textFileArray(:,1) == 270);
% ind315 = find(textFileArray(:,1) == 315);
% 
% 
% 
% %% now tuning curves
% disp('Computing DS'); %Status update
% 
% %Zeros start
% DSI = zeros(num_cells,1);
% vecSum = zeros(num_cells,1);
% vecTheta = zeros(num_cells,1);
% maxDF = zeros(num_cells,1);
% 
% rhos_all = zeros(num_cells, 9);
% rhos1_all = zeros(num_cells, 9);
% rhos2_all = zeros(num_cells, 9);
% rhos3_all = zeros(num_cells, 9);
% rhos_var = zeros(num_cells, 9);
% %Zeroes end
% 

