addpath(genpath('.\GTWgraphConstruct'));
clear all;
clc;

% load demo_data
load('demo_data.mat');

figure;
for i = 1:size(data,3)
    imshow(data(:,:,i),'InitialMagnification','fit');
    pause(.2);
end

%% setting
kappa = 1;
smo = kappa/2;
metric = 'squared';
winSize = T;

%% BILCO
data = reshape(data,[],T);
pix = [1:H*W]';
ref = repmat(ref0,numel(pix),1);
tst = data(pix,:);
[head,tail] = getNeighborRelation(pix,H,W);
Gij = [head,tail];

% get initialCut
distMatrix = zeros(T,T,H*W);
for i = 1:H*W
   distMatrix(:,:,i) = (ref(i,:)-tst(i,:)').^2;
end
initialCut0 = DTW_Edge_input(mean(distMatrix,3))+1;
smo_initial = repmat(initialCut0,H*W,1);
clear distMatrix;

% BILCO
tic;
[midPoints,maxflowBILCO] = BILCO(ref,tst,Gij,smo,smo_initial,winSize,metric);
tBILCO = toc;

% convert to warping paths
paths = cell(H*W,1);
for i = 1:H*W
    paths{i} = midPoint2path(midPoints(i,:),T,T);
end
%% display result
thrVec = 0.4:0.1:0.6;
tDlyCPR = getDelayMidPoints(midPoints,ref,thrVec);
riseMap = nan(H,W);
riseMap(pix) = tDlyCPR;
figure;imagesc(riseMap);    % blue is earlier, yellow is later