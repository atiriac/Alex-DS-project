%% Alex's calculate DSI (and other things) function
% Created 2/7/2020 by Alex Tiriac
%
% Syntax: [DSI,vecSum, vecTheta, wvf_resp_mean, sumVar, rhos_all, rhos_norm_all] = calcDS(roiInt, textFileArray, frameRate)
% ______
% Inputs
% - roiInt: 2D matrix where each row is a cell and each column is a frame's
% fluorescence intensity for that cell.
% - textFileArray: CSV file where each row is a trial. First column needs
% to be direction of moving stim. 7th column needs to be start of stim, 8th
% column needs to be end of stim.
% - frameRate: framerate of movie.
% ______
% Outputs
% - DSI: direction selectivity index of each cell
% - vecSum: vector sum of each cell
% - vecTheta: preferred direction of each cell
% - wvf_resp_mean: chunked mean response of each cell to all directions.
% - sumVar: sum of the variance of each cell, can be used to determine how
% reliable a cell's response to a direction is.
% rhos_all: mean response of each cell at each of the direction



function [DSI,vecSum, vecTheta, wvf_resp_mean, sumVar, rhos_all] = calculateDS(roiInt, textFileArray, frameRate)

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
