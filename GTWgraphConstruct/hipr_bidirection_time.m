function [time] = hipr_bidirection_time(ref,tst,Gij,smo,initial)
%% HIPR + bidirectional strategy
    %% since HIPR didn't provide cut, here use the result of IBFS with window to obtain the cut in one direction
    [N,T] = size(tst);
    [ ss0,ee0,~ ] = build_cGTWGraph( ref, tst, Gij(:,1), Gij(:,2), smo, length(ref), 1);
    
    [mapMatrix] = mapSourceSink( ref, tst, initial, ones(size(initial))*(T-1));
    ss = ss0;
    ss(mapMatrix,1) = 1e8; % label as source
    [maxflowIBFS, label_mid] =  aoIBFS.graphCutMex(ss,ee0);
    paths = label2path(label_mid,N,T);
    midBound = getMidPoints(paths);


    [ ss,ee,~ ] = build_cGTWGraph_low_high( ref, tst, Gij(:,1), Gij(:,2), smo, 1, initial,ones(N,T-1)*(T-1));
    ss = round(ss * 1000);
    ee(:,3:4) = round(ee(:,3:4) * 1000);
    tic;
    maxflowHIPR = hiprMex(ss,ee,N);
    maxflowHIPR = maxflowHIPR/1000;
    time = toc;

    [ ss,ee,~ ] = build_cGTWGraph_low_high( ref, tst, Gij(:,1), Gij(:,2), smo, 1, zeros(N,T-1),midBound);
    ss = ss(:,[2,1]);
    ee = ee(:,[1,2,4,3]);
    ss = round(ss * 1000);
    ee(:,3:4) = round(ee(:,3:4) * 1000);
    tic;
    maxflowHIPR = hiprMex(ss,ee,N);
    maxflowHIPR = maxflowHIPR/1000;
    time2 = toc;
    time = time + time2;
end

