%% findDS
% this code find DS pixels based on a moving 2D avg around every pixel
% 
% Inputs:
% -dF/F movie file
% -DS text file


%% Inputs
frameRate = 1.48;  % 1.48

blurVar = 2;

threshold_DS = 0.13; %0.15 for calDyes   0.2 for g6
threshold_VS = 0.11; %0.11 for calDyes   0.2 for g6
threshold_maxDF = 0;


%Direction correction from DMD comp to rig
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
disp('pick the registered dF/F imaging file'); 
[movie_name, path_name] = uigetfile([path, 'experiment', '\*.tif']);

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
movieBlur = zeros(height, width, n_frames);
for k = 1:n_frames
    movieMat(:,:,k) = imread(movie_name, k);
%     movieBlur(:,:,k) = medfilt2(movieMat(:,:,k),[blurVar, blurVar]); %med
    movieBlur(:,:,k) = imgaussfilt(movieMat(:,:,k),blurVar); %Gauss
end
%End loading the movie

%Apply blur &Comment out if you don't want to blur
movieMat = movieBlur;
%End apply blur

disp('Done loading movie')

%% Correct directions
newInd = zeros(24,1);
replaceDir = [real0; real45; real90; real135; real180; real225; real270; real315];
temp = unique(textFileArray(:,1));

for i = 1:numel(temp)
    indices = textFileArray(:,1) == temp(i);
    newInd(indices) = replaceDir(i);
end
textFileArray(:,1) = newInd;



%% Reshape the movie into one line and compute DS/VS/maxDF
% Use reshape(x, 1, row*col,frames)
% When ready to bring it back into shape, use reshape (y, row, col, frames)
%

movieLine = reshape(movieMat, height*width, n_frames);
[num_cells, dummy] = size(movieLine);

%This is my old DS code
barResp = zeros(num_cells,num_trials);
wvf_resp = zeros(num_cells, 22,num_trials); % waveforms of all responses
wvf_resp_allComb = [];

for i = 1:num_trials
    
    stimFrameStart = floor(textFileArray(i,7)*frameRate);
    stimFrameEnd = ceil(textFileArray(i,8)*frameRate);
    barResp(:,i) = max(movieLine(:,stimFrameStart:stimFrameEnd+3),[],2);
    
    wvf_resp(:,:,i) = movieLine(:,stimFrameStart-3:stimFrameStart+18);
    wvf_resp_allComb = [wvf_resp_allComb; movieLine(:,stimFrameStart-3:stimFrameStart+18)];
    
end

ind0 = find(textFileArray(:,1) == 0);
ind45 = find(textFileArray(:,1) == 45);
ind90 = find(textFileArray(:,1) == 90);
ind135 = find(textFileArray(:,1) == 135);
ind180 = find(textFileArray(:,1) == 180);
ind225 = find(textFileArray(:,1) == 225);
ind270 = find(textFileArray(:,1) == 270);
ind315 = find(textFileArray(:,1) == 315);



%% now tuning curves

%Zeros start
DSI = zeros(num_cells,1);
vecSum = zeros(num_cells,1);
vecTheta = zeros(num_cells,1);
maxDF = zeros(num_cells,1);

rhos_all = zeros(num_cells, 9);
rhos1_all = zeros(num_cells, 9);
rhos2_all = zeros(num_cells, 9);
rhos3_all = zeros(num_cells, 9);
%Zeroes end

% temp_barResp = barResp;
% barResp = barResp - 1; %Uncomment this ifyou want to remove back
% barResp(barResp<0) = 0;
thetas = [0, pi/4, pi/2, 3*pi/4, pi, 5*pi/4, 3*pi/2, 7*pi/4, 0]; %Theta values for 0 45 90 135 180 225 270 315


for i = 1:num_cells
    
    rhos = [mean(barResp(i,ind0)), mean(barResp(i,ind45)), mean(barResp(i,ind90)), mean(barResp(i,ind135)), mean(barResp(i,ind180)), mean(barResp(i,ind225)), mean(barResp(i,ind270)), mean(barResp(i,ind315)), mean(barResp(i,ind0))];
    rhos1 = [barResp(i,ind0(1)), barResp(i,ind45(1)), barResp(i,ind90(1)), barResp(i,ind135(1)), barResp(i,ind180(1)), barResp(i,ind225(1)), barResp(i,ind270(1)), barResp(i,ind315(1)), barResp(i,ind0(1))];
    rhos2 = [barResp(i,ind0(2)), barResp(i,ind45(2)), barResp(i,ind90(2)), barResp(i,ind135(2)), barResp(i,ind180(2)), barResp(i,ind225(2)), barResp(i,ind270(2)), barResp(i,ind315(2)), barResp(i,ind0(2))];
    rhos3 = [barResp(i,ind0(3)), barResp(i,ind45(3)), barResp(i,ind90(3)), barResp(i,ind135(3)), barResp(i,ind180(3)), barResp(i,ind225(3)), barResp(i,ind270(3)), barResp(i,ind315(3)), barResp(i,ind0(3))];
    
  
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


figure('Name', 'cell DS, VS, & max dF', 'position', [100, 100, 1000, 400])
subplot(1,2,1)
plot(vecSum(:,1), DSI(:,1),'k.');
xlabel('Norm Vector Sum')
ylabel('DSI')
subplot(1,2,2)
histogram(maxDF, 'BinLimits', [0, 3], 'BinWidth', 0.05)

%% reshape back into a 2D movie after some thresholds



DSmask = DSI > threshold_DS;
VSmask = vecSum > threshold_VS;
DFmask = maxDF > threshold_maxDF;

mask = DSmask .* VSmask .* DFmask;

thetaImg = reshape(vecTheta, height, width);
maskImg = reshape(mask, height, width);

thetaMasked = thetaImg .* maskImg;

hF = figure('units', 'pixel', 'position', [200, 100, width, height]);
hA1 = axes('Units','Normalized','Position',[0, 0, 1, 1]);
imagesc(thetaMasked)
cmap = jet(1000);
cmap(1,:) = [0,0,0];
colormap(cmap)

saveas(gcf,'whereDSare.jpg')

thetaMaskedFilt = medfilt2(thetaMasked,[3,3]);
hF2 = figure('units', 'pixel', 'position', [200, 100, width, height]);
hA2 = axes('Units','Normalized','Position',[0, 0, 1, 1]);
imagesc(thetaMaskedFilt)
cmap = jet(1000);
cmap(1,:) = [0,0,0];
colormap(cmap)

% J = im2uint8(thetaMasked);
% saveastiff(J, 'test.tif');



% imwrite(thetaMasked,'whereDSare.tif', 'ColorSpace', 'icclab')
% 
% 
% 
% 
% figure, imagesc(thetaMasked)






% 
% %% Light stims
% 
% 
% [sorted, ind] = sort(textFileArray(:,1));
% stimFrameStart = floor(textFileArray(:,7)*frameRate);
% stimFrameEnd = ceil(textFileArray(:,8)*frameRate);
% 
% avgFrame = zeros(height, width, length(ind),1);
% for i = 1:length(ind)
%     avgFrame(:,:,i) = mean(movieMat(:,:,stimFrameStart(ind(i)):stimFrameEnd(ind(i))),3);
% end
% 
%     