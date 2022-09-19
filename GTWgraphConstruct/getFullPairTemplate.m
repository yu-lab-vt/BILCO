function [edge,src,sink] = getFullPairTemplate(T1,T2)
% getFullPairTemplate: construct a DTW dual graph with all nodes
% edge: the first column is startpoint in dual graph
% the second column is end point in dual graph
% the third column is the startpoint in primal graph

    edge = [];
    if(T1==0 || T2==0)
        src = 1;
        sink = 2;
       return; 
    end
    % label of source and sink
    dualPoints = (T1-1)*(T2-1)*2;
    src = dualPoints + 1;
    sink = src + 1;
    labels = zeros(T1-1,T2-1,2);
    labels(:) = 1:dualPoints;
    %% %% three direction edges
    %% right arrow
    % the bottom line
    startPoint = ones(T1-1,1)*src;
    endPoint = [1:T1-1]';
    curWL = [2:T1]';
    edge = [edge;startPoint,endPoint,curWL];
    % upper lines
    startPoint = [dualPoints/2 + 1:dualPoints]';
    endPoint = startPoint - dualPoints/2 + T1-1;
    endPoint(endPoint>dualPoints/2) = sink;
    curWL = repmat([2:T1]',1,T2-1) + T1*repmat(1:T2-1,T1-1,1);
    curWL = curWL(:);
    edge = [edge;startPoint,endPoint,curWL];
    
    %% inclined arrow
    startPoint = [1:dualPoints/2]';
    endPoint = startPoint + dualPoints/2;
    curWL = curWL;
    edge = [edge;startPoint,endPoint,curWL];
    
    %% up arrow
    % the left line
    startPoint = dualPoints/2 + [1:T1-1:dualPoints/2]';
    endPoint = ones(T2-1,1) * sink;
    curWL = [1:T1:T1*(T2-2)+1]'+T1;
    edge = [edge;startPoint,endPoint,curWL];
    % right lines
    endPoint = zeros(T1-1,T2-1);
    endPoint(:) = 1:dualPoints/2;
    startPoint = endPoint + dualPoints/2 + 1;
    startPoint(end,:) = src;
    startPoint = startPoint(:);
    endPoint = endPoint(:);
    curWL = repmat([2:T1]',1,T2-1) + T1*repmat(0:T2-2,T1-1,1);
    curWL = curWL(:)+T1;
    edge = [edge;startPoint,endPoint,curWL];
end