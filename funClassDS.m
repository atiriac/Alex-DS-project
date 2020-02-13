%% k-means functional classification of DS based on preferred directions
% Created by Alex Tiriac on 2/8/2020
% ______
% Inputs
% - prefDir: cell's preferred directions
% - vecSum: cell's vector sum
% - numClusters: number of cluster to classify into, this should be informed
% by a silhouette analysis
%
% ______
% Outputs
% - x: x coordinate of cell's preferred vector
% - y: y coordinate of cell's preferred vector
% - idxOrdered: the functional cluster cell is in 1=T, 2=D, 3=N, 4=V
% - C: the centroid of each functional group




function [x,y,idxOrdered, C] = funClassDS(prefDir,vecSum, numClusters)

    [x,y] = pol2cart(prefDir,vecSum);
    [xNorm,yNorm] = pol2cart(prefDir,ones(length(prefDir),1));

    [idx,C] = kmeans([xNorm yNorm],numClusters);

    idxOrdered = idx;
%     idxOrdered = nan(length(idx),1);
% 
%     % Reorder grps based on max X or Y
%     %Temporal = 1
%     [dummy,id] = max(C(:,1));
%     idxOrdered(idx == id) = 1;
%     %Dorsal = 2
%     [dummy,id] = max(C(:,2));
%     idxOrdered(idx == id) = 2;
%     %Nasal = 3
%     [dummy,id] = min(C(:,1));
%     idxOrdered(idx == id) = 3;
%     %Ventral = 4
%     [dummy,id] = min(C(:,2));
%     idxOrdered(idx == id) = 4;
%     %Done with Reorder

end