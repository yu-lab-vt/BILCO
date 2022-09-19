%% application to extracting depth information from binocular stereo
clear all;
close all;
clc;
addpath(genpath('.\referenced_methods'));
addpath(genpath('.\GTWgraphConstruct'));
p0 = 'D:\experiments data\Binocular stereo\left\';
p1 = 'D:\experiments data\Binocular stereo\right\';
files0 = dir(fullfile(p0,['*.png']));
files1 = dir(fullfile(p1,['*.png']));
rescl = 0.3;
metric = 'absolute';
percent = 0.05;
smos = [0.5,0.5,0.1,0.5,0.5,0.5];

for xx = 1:numel(files0)
    close all;
    %% DTW
    im0 = imresize(double(rgb2gray(imread([p0,files0(xx).name]))),rescl);
    im1 = imresize(double(rgb2gray(imread([p1,files1(xx).name]))),rescl);
    f0 = files0(xx).name;
    f0 = [f0(1:end-4),',mat'];
    smo = smos(xx);
    
    [H,W] = size(im0);
    winSize = round(W/5);
    
    disparityMap = nan(H,W);
    for k = 1:H
       [dist,ix,iy] = dtw(im0(k,:),im1(k,:),metric);
       disparityMap(k,ix) = ix-iy;

       for x = 1:W
           id = find(ix==x);
           if(numel(id)==1)
               disparityMap(k,x) = ix(id) - iy(id);
           elseif(numel(id)>1)
               disparityMap(k,x) = mean(ix(id) - iy(id));
           end   
       end
       preX = 1;
       x = 2;
       while(x<W)
           while(isnan(disparityMap(k,x)))
               x = x + 1;
           end
           disparityMap(k,preX:x) = [1:x-preX+1]*(disparityMap(k,x)-disparityMap(k,preX))/(x-preX) + disparityMap(k,preX);
           preX = x;
           x = x + 1;
       end
    end

    disparityMap(im0==0) = nan;
    %% BILCO
    s = [1:H-1]';
    t = [2:H]';
    Gij = [s,t];
        
    initialCut0 = [1:W-1]+0.5;
    initialCut = repmat(initialCut0,H,1);
    tic;
    [midPoints,maxflowBILCO] = BILCO(im0,im1,Gij,smo,initialCut,winSize,metric);
    tBILCO = toc;

    %% compare disparity map
    disparityMap2 = getDisparityMap(midPoints,H,W);
    disparityMap2(im0==0) = nan;

    figure;
    colormap(gray);
    clims = [quantile(disparityMap2(:),percent),quantile(disparityMap2(:),1-percent)];
    subplot(2,1,1);
    imagesc(disparityMap,clims);
    subplot(2,1,2);
    imagesc(disparityMap2,clims);
    
    %% IBFS
    [ ss,ee,mapMatrix ] = build_cGTWGraph( im0, im1, Gij(:,1), Gij(:,2), smo, winSize, 1, [], metric);
    tic;
    [maxflowIBFS,labelsOrg] = aoIBFS.graphCutMex(ss,ee);
    tIBFS = toc;
    
    %% BK
    tic;
    maxflowBK = graphCutDynamicMex(ss, ee);
    tBK = toc;
    
    %% HIPR
    ss = round(ss * 10);
    ee(:,3:4) = round(ee(:,3:4) * 10);
    %% HIPR
    tic;
    maxflowHIPR = hiprMex(ss,ee,H*W);
    maxflowHIPR = maxflowHIPR/10;
    tHIPR = toc;
    
    %% HPF
    nNode = size(ss,1);
    src = nNode+1;
    sink = nNode+2;
    col3 = [ss(:,1);ss(:,2);ee(:,3);ee(:,4)];
    clear ss;
    ee = ee(:,1:2);
    col1 = [ones(nNode,1)*src;[1:nNode]';ee(:,1);ee(:,2)];
    col2 = [[1:nNode]';ones(nNode,1)*sink;ee(:,2);ee(:,1)];
    clear ee;
    mat = sparse(col1,col2,col3,nNode+2,nNode+2);

    clear col1;
    clear col2;
    clear col3;
    tic;
    maxflowHPF = hpf(mat,src,sink);
    maxflowHPF = maxflowHPF/10;
    tHPF = toc;
end