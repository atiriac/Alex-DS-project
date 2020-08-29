%% Build the neuron table
% this code will build a table where every row is a neuron
% Created by Alex T on 3/13/2019
%
% This code needs an excel spreadsheet with metadata in it.
% Right now, it will read 'DSmainExperimentsSheet.xlsx'
% Every row in that excell spreadsheet is a field of view.
% The buildNeuronTable code will run all the field of views
%
% The end result of this code is that it saves 'neuronTable.mat'
% This is the table where each row is a neuron
% Other codes will load this table and run analysis on it

%% Load the table
data_guide_name = 'DSmainExperimentsSheet.xlsx';

% Detect import options for the data guide spreadsheet
opts = detectImportOptions(data_guide_name);

% Ensure that the 'data_name' column of the table is a string array
opts = setvartype(opts,{'ExperimentDate', 'FileName','LightCondition','GFPLabel_','eye', 'location', 'On_OffCells',...
    'OnCells', 'BadCells','GFPcellID','CalciumSensor'},'string');

expTable = readtable(data_guide_name, opts');

[num_files dummy] = size(expTable);

%% Intialize the neuron table (called neuronTable)

totalNeurons = sum(expTable.numberOfNeurons);
neuronCounter = 1; %This value will keep track of which neuron I am on

neuronTable = table('Size', [totalNeurons 24], 'VariableTypes', {'double','double','string','string','string',...
    'double','string','string','string','string','string','double','double','double','string','double','double',...
    'string','double','double','double','double','double','double'});

neuronTable.Properties.VariableNames = {'neuronNum', 'neuronNumWithinFOV', 'animalID','fileName', 'sex',...
    'age', 'condition', 'calciumSensor','GFPretina','GFPid','eye','corrX','corrY','distanceFromON','location','degCorr','frameRate',...
    'cellID','DSI','DSIsig','vecSum','prefDir','prefDirCorr','varSum'};

neuronTable.wvfRespToBars = NaN*ones([totalNeurons,160]);
neuronTable.meanRespToBars = NaN*ones([totalNeurons,9]);
neuronTable.meanRespToBarsNormalized = NaN*ones([totalNeurons,9]);
neuronTable.cellLoc = NaN*ones([totalNeurons,2]);
neuronTable.shuffledDSI = NaN*ones([totalNeurons,1000]);

neuronTable.GFPid = repmat("none",[totalNeurons 1]);


%% Here is the for loop that goes through every file


for i = 1:num_files   % this will run through all files
% for i = 1:2 % For testing code

    %load file(i) and initialize some variables
    load(expTable.FileName(i)) 
    animalID = char(expTable.FileName(i));
    animalID = animalID(1:6);
    
    % Call function calcDS to calculate DS
    [DSI, vecSum, vecTheta, wvf_resp_mean, sumVar, rhos_all, rhos_norm_all] = calcDS(roiInt, textFileArray, expTable.FrameRate(i));
   
    % Determine if DSI is significant via bootstrapping
    shuffledDSI = nan(length(DSI),1000);
    for j = 1:1000
        shuffledText = textFileArray;
        shuffledText(:,1) = shuffleTrialDirections(textFileArray(:,1),8, 3);

        [shuffledDSI(:,j), dummy, dummy, dummy, dummy, dummy, dummy] = calcDS(roiInt, shuffledText, expTable.FrameRate(i));
    end
    
    % Call function to make the CellID string array
    [cellID, numID] = determineID(str2num(expTable.On_OffCells(i)), str2num(expTable.OnCells(i)), str2num(expTable.BadCells(i)),expTable.numberOfNeurons(i));
    
    % Call function to add labels to GFPid so we can identify GFP+ cells
    if expTable.GFPLabel_(i) ~= "none"
        [GFPid] = determineGFPid(str2num(expTable.GFPcellID(i)),expTable.numberOfNeurons(i),expTable.GFPLabel_(i));
        neuronTable.GFPid(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = GFPid;
    end
    

    % Load neuron info in the neuronTable
    neuronTable.neuronNum(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = [neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1]';
    neuronTable.neuronNumWithinFOV(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = [1:expTable.numberOfNeurons(i)]';
    neuronTable.animalID(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(animalID,[expTable.numberOfNeurons(i),1]);
    neuronTable.fileName(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.FileName(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.sex(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.Sex(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.age(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.Age(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.condition(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.LightCondition(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.calciumSensor(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.CalciumSensor(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.GFPretina(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.GFPLabel_(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.eye(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.eye(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.corrX(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.Corr_x(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.corrY(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.Corr_y(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.distanceFromON(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.DistanceFromON_um_(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.location(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.location(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.degCorr(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.degCorrection(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.frameRate(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = repmat(expTable.FrameRate(i),[expTable.numberOfNeurons(i),1]);
    neuronTable.cellID(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = cellID;
    neuronTable.DSI(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = DSI;
    neuronTable.shuffledDSI(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = shuffledDSI;
    neuronTable.vecSum(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = vecSum;
    neuronTable.prefDir(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = vecTheta;
    neuronTable.varSum(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1) = sumVar;
    neuronTable.wvfRespToBars(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = wvf_resp_mean;
    neuronTable.meanRespToBars(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = rhos_all;
    neuronTable.meanRespToBarsNormalized(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = rhos_norm_all;

    if exist('cellCenters','var') == 1
        neuronTable.cellLoc(neuronCounter:neuronCounter+expTable.numberOfNeurons(i)-1,:) = cellCenters;
        clear cellCenters
    end
    
    %Update how many neurons we have analyzed so far (+1)
    neuronCounter = neuronCounter+expTable.numberOfNeurons(i);
    
    if floor(i/10)==i/10
        i
    end
end

neuronTable.prefDirCorr = neuronTable.prefDir + deg2rad(neuronTable.degCorr);
for i = 1:totalNeurons
    %Fix prefDirCorr
    if neuronTable.prefDirCorr(i) > 2*pi
        neuronTable.prefDirCorr(i) = neuronTable.prefDirCorr(i) - 2*pi;
    elseif neuronTable.prefDirCorr(i) < 0
        neuronTable.prefDirCorr(i) = neuronTable.prefDirCorr(i) + 2*pi;
    end
    
    %Determine DSI significance
    [val,idx]=min(abs(neuronTable.DSI(i)-sort(neuronTable.shuffledDSI(i,:))));
    neuronTable.DSIsig(i) = idx/1000;
end

%classify all non-BAD hb9 and drd4 cells as ON-OFF cells
neuronTable.cellID(strcmp(neuronTable.GFPid, 'hb9')==1 & strcmp(neuronTable.cellID,'BAD')==0) = 'ON-OFF';
neuronTable.cellID(strcmp(neuronTable.GFPid, 'drd4')==1 & strcmp(neuronTable.cellID,'BAD')==0) = 'ON-OFF';


save('neuronTable.mat','neuronTable');









%% Functions

%DS function
function [DSI,vecSum, vecTheta, wvf_resp_mean, sumVar, rhos_all, rhos_norm_all] = calcDS(roiInt, textFileArray, frameRate)

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


% Build a string identifier for cells (Are they ON or ON-OFF?)
function [cellID, numID] = determineID(onOffCells, onCells, badCells, numNeurons)
    
numID = zeros(numNeurons,1);
numID(onOffCells) = 1;
numID(onCells) = 2;
for i = 1:numNeurons
    if numID(i) == 1
        cellID(i,1) = "ON-OFF";
    elseif numID(i) == 2
        cellID(i,1) = "ON";
    else
        cellID(i,1) = "BAD";
    end
end

end


function [GFPid] = determineGFPid(indexGFP,numNeurons,GFPname)

numID = zeros(numNeurons,1);
numID(indexGFP) = 1;
for i = 1:numNeurons
    if numID(i) == 1
        GFPid(i,1) = GFPname;
    else
        GFPid(i,1) = "none";
    end
end

end