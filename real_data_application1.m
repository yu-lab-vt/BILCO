%% application to calculating signal propagation
close all;
clear all;
clc;
addpath(genpath('.\referenced_methods'));
addpath(genpath('.\GTWgraphConstruct'));
%% setting
p0 = 'D:\experiments data\Signal propagation\';
files = dir(fullfile(p0,['*.mat']));

for ii = 1:numel(files)
    f0 = files(ii).name;
    load([p0,f0]);
    ref = double(ref);
    tst = double(tst);
    
    %% BILCO
    T = size(tst,2);
    distMatrix = zeros(T,T,numel(pix));
    for i = 1:numel(pix)
       distMatrix(:,:,i) = (ref(i,:)-tst(i,:)').^2;
    end
    initialCut0 = DTW_Edge_input(mean(distMatrix,3))+1;
    initialCut = repmat(initialCut0,numel(pix),1);
    clear distMatrix;

    tic;
    [midPoints,maxflowBILCO] = BILCO(ref,tst,Gij,smo,initialCut);
    tBILCO = toc;
    
    %% IBFS
    [ ss,ee,mapMatrix ] = build_cGTWGraph( ref, tst, Gij(:,1), Gij(:,2), smo, length(ref), 1);
    tic;
    [maxFlowIBFS,labelsOrg] = aoIBFS.graphCutMex(ss,ee);
    tIBFS = toc;
    
    %% BK
    tic;
    maxflowBK = graphCutDynamicMex(ss, ee);
    tBK = toc;
    
    ss = round(ss * 1000);
    ee(:,3:4) = round(ee(:,3:4) * 1000);
    %% HIPR
    tic;
    maxflowHIPR = hiprMex_large_dataset(ss,ee,numel(pix));
    maxflowHIPR = maxflowHIPR/1000;
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
    maxflowHPF = maxflowHPF/1000;
    tHPF = toc;
end