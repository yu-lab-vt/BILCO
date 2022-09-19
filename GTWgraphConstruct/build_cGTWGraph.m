function [ ss,ee,mapMatrix ] = build_cGTWGraph( ref, tst, s, t, smoBases, maxShift, stds, path0, metric)
% build_cGTWGraph: cGTW graph construction by allowing shift on estimated
% path
% Use coordinate to specify nodes in template pair, but use integer for spatial
% Support arbitrary graph

% INPUT
% ref: 1xT or NxT reference curve
% tst: NxT curves for pixels
% s,t: graph edge pair, undirected, if has pair (i,j), do not include (j,i), do not use (i,i) 
% smoBases: smoothness between edges
% maxShift: window size, (winSize-1) lines above and below diagonal
% stds: noise variance for each test curve, or use single s2 for all test curves
% path0: estimated warping path. If no path0, it is selected as the diagonal.

% OUTPUT
% ss:  -	the edges connecting the source and the sink with the regular nodes (array of type double, size : [numNodes, 2])
% 				termWeights(i, 1) is the weight of the edge connecting the source with node #i
% 				termWeights(i, 2) is the weight of the edge connecting node #i with the sink
% 				numNodes is determined from the size of termWeights.
% ee:  - 	the edges connecting regular nodes with each other (array of type double, array size [numEdges, 4])
% 				edgeWeights(i, 3) connects node #edgeWeights(i, 1) to node #edgeWeights(i, 2)
% 				edgeWeights(i, 4) connects node #edgeWeights(i, 2) to node #edgeWeights(i, 1)
%				The only requirement on edge weights is submodularity: edgeWeights(i, 3) + edgeWeights(i, 4) >= 0

% Basic constraint:
% No skipping/turn left/go down: otherwise has infinite cost (capacity)
% Start NEAR (1,1) and stop NEAR (T,T)
%
% The graphs has x axis as REF and y axis as TST
%
% ALL INPUTS SHOULD BE DOUBLE

T1 = size(ref,2);
[nNode,T2] = size(tst);
cap = 1e8;


if(~exist('path0','var') || isempty(path0))
    path0 = cell(nNode,1);
    for i = 1:nNode
        if(T1<=T2)
            path0{i}  = [repmat([1:T1]',1,2);ones(T2-T1,1)*T1,[T1+1:T2]'];
        else
            path0{i}  = [repmat([1:T2]',1,2);,[T2+1:T1]',ones(T1-T2,1)*T2];
        end
    end
end

if(~exist('metric','var'))
    metric = 'squared';
end

if(numel(smoBases)==1)
   smoBases = ones(nNode,1) *smoBases;
end
if(numel(stds)==1)
    stds = ones(nNode,1) *stds;
end

if(size(ref,1)==1)
    ref = repmat(ref,nNode,1);
end

%%  Full DTW dual graph template
% return the edges in dual graph, and start point of each edge in primal
% graph for mapping the weight
[TemplateEdge,source,sink] = getFullPairTemplate(T1,T2);
% set each DTW graph has one substitute source and sink.
% first, set each DTW graph and its dual graph
% still, use coordinates

% the matrix which record the mapping (Or the status whether the node is valid)
nNodeInDualDTW = (T1-1)*(T2-1)*2+2;
mapMatrix = zeros(nNode,nNodeInDualDTW);

% tic;
%% the dual node coordinate relative in primal graph
[x,y,z] = ind2sub([T1-1,T2-1,2],(1:(T1-1)*(T2-1)*2)');
y_location = y - 0.5 + (z-1.5)/2;
nEdgeInGraph = size(TemplateEdge,1);
ss = zeros(nNode*nNodeInDualDTW,2);
ee = zeros(nNode*nEdgeInGraph,4);
nValidEdges = 0;
for k = 1:nNode
    % allow shift on the estimated warping path
   [lBcurve,uBcurve] = dilatePath(path0{k},maxShift);
   mapToSrcSink = 1:nNodeInDualDTW;
   %% set validDualGraph
   mapToSrcSink(y_location<lBcurve(x)') = source;
   mapToSrcSink(y_location>uBcurve(x)') = sink;
   
   %% %%  Plan for faster construction
   % 1. for each DTW dual graph, construct the full graph, that is one template
   % 2. set non-valid dual nodes to source or sink
   % 3. remove the edges which both start and end are source or sink
   %% replace the non-valid nodes to source or sink.
   curEdge = TemplateEdge;
   curEdge(:,1) = mapToSrcSink(curEdge(:,1));
   curEdge(:,2) = mapToSrcSink(curEdge(:,2));
   %% remove non-valid edges (both points are source or sink)
   validEdges = curEdge(:,2)-curEdge(:,1)~=0;
   curEdge = curEdge(validEdges,:);
   %% weights
   d0 = (getDistMat(ref(k,:),tst(k,:),metric)/stds(k));
   weights = d0(curEdge(:,3));
   %% Update
   % link substitute source and sink to primal source and sink
   ss(source + (k-1)*nNodeInDualDTW,:) = [cap,0];
   ss(sink + (k-1)*nNodeInDualDTW,:) = [0,cap];
   % update node label
   curEdge(:,1:2) = curEdge(:,1:2) + (k-1)*nNodeInDualDTW;
   % update ee
   nCurValidEdge = size(curEdge,1);
   ee(nValidEdges + 1:nCurValidEdge + nValidEdges,:) = [curEdge(:,1:2),weights,ones(numel(weights),1)*cap];
   nValidEdges = nValidEdges + nCurValidEdge;
   % update others
   mapMatrix(k,:) = mapToSrcSink;
end
ee = ee(1:nValidEdges,:);
% toc;
% tic;
if(sum(smoBases)>0)
    %% neighbor dual graph Linking.
    % for each neighbor pair, set bidirected edges with smoBase weight
    eeSpa = zeros(nNodeInDualDTW*numel(s),4);
    validEdges = false(nNodeInDualDTW*numel(s),1);
    for k = 1:numel(s)
        id1 = s(k);
        id2 = t(k);
        smoBase = max(smoBases(id1),smoBases(id2));
        curEE = [mapMatrix(id1,:)',mapMatrix(id2,:)',ones(nNodeInDualDTW,1)*smoBase,ones(nNodeInDualDTW,1)*smoBase];
        indexRange = (k-1)*nNodeInDualDTW+1 : k*nNodeInDualDTW;
        validEdges(indexRange) = curEE(:,1)<source | curEE(:,2)<source | curEE(:,1)~=curEE(:,2);
        curEE(:,1) = curEE(:,1) + (id1-1)*nNodeInDualDTW;
        curEE(:,2) = curEE(:,2) + (id2-1)*nNodeInDualDTW;
        eeSpa(indexRange,:) = curEE;
    end
    eeSpa = eeSpa(validEdges,:);
    ee = [ee;eeSpa];
end

nTotalNode = numel(ss);
[uniqueNodes] = unique(ee(:,1:2));
validNodes = false(nTotalNode,1); 
validNodes(uniqueNodes) = true;
ss = ss(validNodes,:);
% change edge node labels
mapNodes = 1:numel(ss);
mapNodes(uniqueNodes) = 1:numel(uniqueNodes);
ee(:,1) = mapNodes(ee(:,1));
ee(:,2) = mapNodes(ee(:,2));

% toc;
end