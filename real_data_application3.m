%% application to signature identification
clear all;
close all;
clc;
addpath(genpath('.\referenced_methods'));
addpath(genpath('.\GTWgraphConstruct'));
%% setting
p0 = 'D:\experiments data\Signature online data\MCYT_Signature_100\';
checkType = 'v';
outName = '_result.mat';
N = 5;
metric = 'squared';
Gij = [];
for i = 1:N
   for j = i+1:N
       Gij = [Gij;i,j];
   end
end

tBILCO = zeros(100,25);
tIBFS = zeros(100,25);
tBK = zeros(100,25);
tHIPR = zeros(100,25);
tHPF = zeros(100,25);

for xx = 1:20
    folderID = xx-1;
    folder = num2str(folderID,'%04d');
    filename = [p0,folder,'\',folder,'v',num2str(0,'%02d'),'.fpg'];
    [x,y,z,az,in,pps]=FPG_Signature_Read(filename, 0, 0);
    ref = [normalize(x),normalize(y),normalize(z),normalize(az),normalize(in)]';
    T1 = numel(x);
    
    for i = 1:25
        fileID = i-1;
        filename = [p0,folder,'\',folder,checkType,num2str(fileID,'%02d'),'.fpg'];
        [x,y,z,az,in,pps]=FPG_Signature_Read(filename, 0, 0);
        tst = [normalize(x),normalize(y),normalize(z),normalize(az),normalize(in)]';
        T2 = numel(x);

        distMatrix = zeros(T2,T1,N);
        for k = 1:N
           distMatrix(:,:,k) = (ref(k,:)-tst(k,:)').^2;
        end
        initialCut0 = DTW_Edge_input(mean(distMatrix,3))+1;
        initialCut = repmat(initialCut0,N,1);
        smo = 0.1;
        
        %% DTW
        [midPoints0,maxflowTest] = BILCO(ref,tst,Gij,0,[],max(T1,T2),metric);
        
        %% BILCO
        tic;
        [midPoints,maxflowBILCO] = BILCO(ref,tst,Gij,smo,initialCut,max(T1,T2),metric);
        tBILCO(xx,i) = toc;

        %% IBFS
        [ ss,ee,mapMatrix ] = build_cGTWGraph( ref, tst, Gij(:,1), Gij(:,2), smo, max(T1,T2), 1, [], metric);
        tic;
        [maxFlowIBFS,labelsOrg] = aoIBFS.graphCutMex(ss,ee);
        tIBFS(xx,i) = toc;

        %% BK
        tic;
        maxflowBK = graphCutDynamicMex(ss, ee);
        tBK(xx,i) = toc;

        ss = round(ss * 10);
        ee(:,3:4) = round(ee(:,3:4) * 10);
        %% HIPR
        tic;
        maxflowHIPR = hiprMex(ss,ee,N);
        maxflowHIPR = maxflowHIPR/10;
        tHIPR(xx,i) = toc;

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
        tHPF(xx,i) = toc;
    end

end


function x = normalize(x)
    x = x - mean(x);
    x = x/std(x);
end