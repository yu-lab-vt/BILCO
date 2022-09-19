function [mapMatrix] = mapSourceSink( ref, tst, lB, uB)
T1 = size(ref,2);
[nNode,T2] = size(tst);
cap = 1e8;

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
for k = 1:nNode
   % allow shift on the estimated warping path
   lBcurve = lB(k,:);
   uBcurve = uB(k,:);
   mapToSrcSink = 1:nNodeInDualDTW;
   %% set validDualGraph
   mapToSrcSink(y_location<lBcurve(x)') = source;
   mapToSrcSink(y_location>uBcurve(x)') = sink;
   mapMatrix(k,:) = mapToSrcSink;
end
% toc;
mapMatrix = mapMatrix';
mapMatrix = mapMatrix(:);
mapMatrix = mapMatrix==source;
end