%% synthetic data, mimic signal propagation
addpath(genpath('.\referenced_methods'));
addpath(genpath('.\GTWgraphConstruct'));
clear all;
clc;

%% setting
T = 20;
H = 20;
W = H;
kappa = 0.2;
smo = kappa/2;
SNR = 10;
metric = 'squared';
T1 = T;
T2 = T;
winSize = T;

for kk = 1:100
    
    %% signal propagation from center to boundary
    centerT = round(T/2);
    centerX = H/2;
    centerY = W/2;
    sigma = T/8;
    propspeedX = T/H/sqrt(2);
    propspeedY = T/H/sqrt(2);
    
    % noise-free data generation
    syntheticdata_no_noise = zeros(H,W,T);
    ref0 = zeros(1,T);
    for i = 1:T     % Gaussian shape curve
        ref0(i) = exp(-(i-centerT)^2/2/sigma^2);
    end
    TW = find(ref0>0.01);   % select the part larger than 1% for calculationg noise level
    n_std = sqrt(mean(ref0(TW).^2)/10^(SNR/10));
    ref0 = ref0 / n_std;
    base = -(sqrt(((W-centerY)*propspeedY)^2 + ((H-centerX)*propspeedX)^2))/2;
    for y = 1:W
        for x = 1:H
            shift = base + sqrt(((y-centerY)*propspeedY)^2 + ((x-centerX)*propspeedX)^2);   % shift of current position
            shift = round(shift);
            if(shift<0)
               curve = [ref0(-shift+1:end),zeros(1,min(T,-shift))];
            else
               curve = [zeros(1,min(T,shift)),ref0(1:end-shift)];
            end
            syntheticdata_no_noise(x,y,:) = curve;
        end
    end

    load('random_Seed.mat');
    rng(s);
    
    data = syntheticdata_no_noise + randn(size(syntheticdata_no_noise));
    data = reshape(data,[],T);
    pix = [1:H*W]';
    ref = repmat(ref0,numel(pix),1);
    tst = data(pix,:);
    [head,tail] = getNeighborRelation(pix,H,W);
    Gij = [head,tail];
    
    %% get initialCut
    distMatrix = zeros(T2,T1,H*W);
    for i = 1:H*W
       distMatrix(:,:,i) = (ref(i,:)-tst(i,:)').^2;
    end
    initialCut0 = DTW_Edge_input(mean(distMatrix,3))+1;
    smo_initial = repmat(initialCut0,H*W,1);
    clear distMatrix;
    
    %% BILCO
    tic;
    [midPoints,maxflowBILCO] = BILCO(ref,tst,Gij,smo,smo_initial,winSize,metric);
    tBILCO = toc;
    
    %% ELCO
    initialCut = ones(H*W,T1-1);
    tic;
    [midPoints,maxflowELCO] = BILCO(ref,tst,Gij,smo,initialCut,winSize,metric);
    tELCO = toc;
    
    %% generate GTW graph
    [ ss,ee,~ ] = build_cGTWGraph( ref, tst, Gij(:,1), Gij(:,2), smo, winSize, 1,[],metric);
    
    %% IBFS
    tic;
    [maxflowIBFS,labelsOrg] = aoIBFS.graphCutMex(ss,ee);
    tIBFS = toc;
    
    %% BK
    tic;
    maxflowBK = graphCutDynamicMex(ss, ee);
    tBK = toc;
    
    %% HIPR and HPF could only deal with integer
    ss = ss * 1000;
    ee(:,3:4) = ee(:,3:4) * 1000;
    
    %% HIPR
    tic;
    maxflowHIPR = hiprMex(ss,ee,H*W);
    maxflowHIPR = maxflowHIPR/1000;
    tHIPR = toc;
    
    %% BI-HIPR
    [tBIHIPR] = hipr_bidirection_time(ref,tst,Gij,smo,smo_initial);

    
    %% HPF
    nNode = size(ss,1);
    src = nNode+1;
    sink = nNode+2;
    col1 = [ones(nNode,1)*src;[1:nNode]';ee(:,1);ee(:,2)];
    col2 = [[1:nNode]';ones(nNode,1)*sink;ee(:,2);ee(:,1)];
    col3 = [ss(:,1);ss(:,2);ee(:,3);ee(:,4)];
    mat = sparse(col1,col2,col3,nNode+2,nNode+2);
    
    tic;
    maxflowHPF = hpf(mat,src,sink);
    maxflowHPF = maxflowHPF/1000;
    tHPF = toc;
    
    
    paths = cell(H*W,1);
    for k = 1:H*W
        paths{k} = midPoint2path(midPoints(k,:),T1,T2);
    end
    
    
    %% compare
    if(abs(maxflowBILCO-maxflowIBFS)/maxflowIBFS>1e-3)
        disp('error');
        keyboard;
    else
        disp(['correct ',num2str(kk),' in ',num2str(1000)]);
    end
end
disp('All correct!');