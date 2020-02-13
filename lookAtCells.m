%% Calcium imaging GUI for looking at cells from the singularity mapping exp


%You will need:
%1) your DS calcium imaging movie (needs to be 8 direction repeated 3
%times)
%2) the on-off segregation movie (for example, I do long bars in 4
%directions, no repeats). (needs to be same width and height as 1)
%3) A binary mask of cells of interest (needs to be same width and height
%as calcium movies)
%4) text file with directions and times of stim for movie 1
%5) text file for movie 2



%% Inputs

% Calibrate the directions on the text file
real0 = 270;
real45 = 225;
real90 = 180;
real135 = 135;
real180 = 90;
real225 = 45;
real270 = 0;
real315 = 315;

frameRate = 1.48;  %1.48




%% Compare DS pre and post Training

%Loading the calcium movies
disp('pick the DS calcium imaging file'); 
[movie_DS, path_name] = uigetfile([path, 'experiment', '\*.tif']);

disp('pick the ON-OFF calcium imaging file'); 
[movie_onOff, path_name] = uigetfile([path, 'experiment', '\*.tif']);

%Loading the mask
disp('pick the calcium ROI mask file'); 
[roi_mask_file, path_name] = uigetfile([path, 'experiment', '\*.tif']);

%Loading the text files
disp('pick the text file with DS directions')
[textF_DS, path_name_text] = uigetfile([path, 'experiment', '\*.txt']);

disp('pick the text file for the ON-OFF segregation movie')
[textF_onOff, path_name_text] = uigetfile([path, 'experiment', '\*.txt']);

textFile_DS = readtable(textF_DS);
textFile_onOff = readtable(textF_onOff);

textFileArray_DS = table2array(textFile_DS);
textFileArray_onOff = table2array(textFile_onOff);

clear textF_DS textFile_DS textF_onOff textFile_onOff path_name path_name_text

%Loading the movies
nFrames = numel(imfinfo(movie_DS));
nFramesOnOff = numel(imfinfo(movie_onOff));
[height, width] = size(imread(movie_DS,1));
movieDS = zeros(height, width, nFrames);
movieOnOff = zeros(height, width, nFramesOnOff);
for k = 1:nFrames
    movieDS(:,:,k) = imread(movie_DS, k);
end
for k = 1:nFramesOnOff
    movieOnOff(:,:,k) = imread(movie_onOff, k);
end
%End loading the movies

clear movie_DS movie_onOff

%Loading the mask
roi_mask = imread(roi_mask_file);
countmask=bwlabel(roi_mask);
%End loading the mask



%% Correct directions
newInd = zeros(24,1);
replaceDir = [real0; real45; real90; real135; real180; real225; real270; real315];
temp = unique(textFileArray_DS(:,1));

for i = 1:numel(temp)
    indices = textFileArray_DS(:,1) == temp(i);
    newInd(indices) = replaceDir(i);
end
textFileArrayOld = textFileArray_DS;
textFileArray_DS(:,1) = newInd; 

clear real0 real45 real90 real135 real180 real225 real270 real315


%% use ROI mask to find z profile at each cell

for i=1:nFrames
    partpropsDS=regionprops(countmask,movieDS(:,:,i),'Area','MeanIntensity','Centroid'); 
    for j = 1:max(max(countmask))
        roiIntDS(j,i) = partpropsDS(j).MeanIntensity;
    end
end

for i=1:nFramesOnOff
    partpropsOnOff=regionprops(countmask,movieOnOff(:,:,i),'Area','MeanIntensity','Centroid'); 
    for j = 1:max(max(countmask))
        roiIntOnOff(j,i) = partpropsOnOff(j).MeanIntensity;
    end
end


[nCells, dummy] = size(roiIntDS);

%% calcDS

[DSI,vecSum, vecTheta,wvf_resp_reordered, wvf_resp_mean, sumVar, rhos_all,rhos1_all,rhos2_all,rhos3_all, rhos_norm_all,thetas] = calcDS(roiIntDS, textFileArray_DS, frameRate);


%% Interactive plot


wvf_resp_t1 = []; %These 6 lines just creates the variables
wvf_resp_t2 = [];
wvf_resp_t3 = [];

for i = 1:8 %This for loop will create wvf for each trial (3 trials here)
    wvf_resp_t1 = [wvf_resp_t1, wvf_resp_reordered(:,:,i*3-2)];
    wvf_resp_t2 = [wvf_resp_t2, wvf_resp_reordered(:,:,i*3-1)];
    wvf_resp_t3 = [wvf_resp_t3, wvf_resp_reordered(:,:,i*3)];
end


% Figures to scroll through cells


cellNumb = 1;

% Create figure windows and save the handles 
fig_info = figure; 
set(gcf, 'Position', [10   100   801   100]);
fig_dfof = figure; 
set(gcf, 'Position', [10   200   801   200]);
fig_wvfs = figure; 
set(gcf, 'Position', [10   430   600   200]);
fig_tuning = figure; 
set(gcf, 'Position', [611   430   200   200]);
fig_onOffRawDf = figure; 
set(gcf, 'Position', [10   715   401   200]);
fig_onOffAvgDf = figure; 
set(gcf, 'Position', [411   715   401   200]);


% Create wrapper function for the redraw function called by video_fig so
% that I can pass in 'constant' arguments in addition to updating frame
% argument
lookAtCells_redraw_wrapper = @(cellNumb) lookAtCells_redraw(cellNumb, fig_info,...
    fig_dfof, fig_wvfs,fig_tuning,fig_onOffRawDf,fig_onOffAvgDf,thetas, ...
    roiIntDS,roiIntOnOff,DSI,vecSum, vecTheta,wvf_resp_reordered, wvf_resp_mean, sumVar, ...
    rhos_all,rhos1_all,rhos2_all,rhos3_all, rhos_norm_all,...
    wvf_resp_t1,wvf_resp_t2,wvf_resp_t3); 

% Initialize to the start of the movie
lookAtCells_redraw(1, fig_info,...
    fig_dfof, fig_wvfs,fig_tuning,fig_onOffRawDf,fig_onOffAvgDf,thetas, ...
    roiIntDS,roiIntOnOff,DSI,vecSum, vecTheta,wvf_resp_reordered, wvf_resp_mean, sumVar, ...
    rhos_all,rhos1_all,rhos2_all,rhos3_all, rhos_norm_all,...
    wvf_resp_t1,wvf_resp_t2,wvf_resp_t3); 

% Hand control over to the user
[scroller, axes_handle, scroll_bar_handles, scroll_func] = videofig(nCells, lookAtCells_redraw_wrapper);

set(scroller,'units','normalized','position',[0.4,0.8,0.2,0.1])
























%% functions

%DS function
function [DSI,vecSum, vecTheta,wvf_resp_reordered, wvf_resp_mean, sumVar, rhos_all,rhos1_all,rhos2_all,rhos3_all, rhos_norm_all,thetas] = calcDS(roiInt, textFileArray, frameRate)

    [num_trials dummy] = size(textFileArray);
    [num_cells num_frames] = size(roiInt);

    % Find max intensity at every bar presentation (not discriminating between ON and OFF)

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

    %Zeros start
    DSI = zeros(num_cells,1);
    vecSum = zeros(num_cells,1);
    vecTheta = zeros(num_cells,1);
    maxDF = zeros(num_cells,1);
    sumVar = zeros(num_cells,1);

    rhos_all = zeros(num_cells, 9);
    rhos1_all = zeros(num_cells, 9);
    rhos2_all = zeros(num_cells, 9);
    rhos3_all = zeros(num_cells, 9);
    rhos_var = zeros(num_cells, 9);
    rhos_norm_all = zeros(num_cells, 9);
    %Zeroes end

    thetas = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 0]; %Theta values for 0 45 90 135 180 225 270 315

    for i = 1:num_cells

        rhos = [mean(barResp(i,ind0)), mean(barResp(i,ind45)), mean(barResp(i,ind90)), mean(barResp(i,ind135)), mean(barResp(i,ind180)), mean(barResp(i,ind225)), mean(barResp(i,ind270)), mean(barResp(i,ind315)), mean(barResp(i,ind0))];
        rhos1 = [barResp(i,ind0(1)), barResp(i,ind45(1)), barResp(i,ind90(1)), barResp(i,ind135(1)), barResp(i,ind180(1)), barResp(i,ind225(1)), barResp(i,ind270(1)), barResp(i,ind315(1)), barResp(i,ind0(1))];
        rhos2 = [barResp(i,ind0(2)), barResp(i,ind45(2)), barResp(i,ind90(2)), barResp(i,ind135(2)), barResp(i,ind180(2)), barResp(i,ind225(2)), barResp(i,ind270(2)), barResp(i,ind315(2)), barResp(i,ind0(2))];
        rhos3 = [barResp(i,ind0(3)), barResp(i,ind45(3)), barResp(i,ind90(3)), barResp(i,ind135(3)), barResp(i,ind180(3)), barResp(i,ind225(3)), barResp(i,ind270(3)), barResp(i,ind315(3)), barResp(i,ind0(3))];
        rhos_var(i,:) = [var(barResp(i,ind0)), var(barResp(i,ind45)), var(barResp(i,ind90)), var(barResp(i,ind135)), var(barResp(i,ind180)), var(barResp(i,ind225)), var(barResp(i,ind270)), var(barResp(i,ind315)), var(barResp(i,ind0))];
        sumVar(i) = sum(rhos_var(i,:));
        

        
        
        
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
        rhos_norm_all(i,:) = rhosNorm;
        %End store for cells


    end
    
    
    reordered_ind = [ind0 ; ind45 ; ind90 ; ind135 ; ind180 ; ind225 ; ind270 ; ind315];
    wvf_resp_reordered = wvf_resp(:,:,reordered_ind);
    wvf_resp_mean = [mean(wvf_resp_reordered(:,:,1:3),3),mean(wvf_resp_reordered(:,:,4:6),3),mean(wvf_resp_reordered(:,:,7:9),3),mean(wvf_resp_reordered(:,:,10:12),3),mean(wvf_resp_reordered(:,:,13:15),3),mean(wvf_resp_reordered(:,:,16:18),3),mean(wvf_resp_reordered(:,:,19:21),3),mean(wvf_resp_reordered(:,:,22:24),3)];
    
end





