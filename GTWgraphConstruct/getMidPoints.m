function midPoints = getMidPoints(paths)
% getMidPoints: get the middle height for each column

    N = numel(paths);
    T = paths{1}(end,1);
    midPoints = zeros(N,T-1);
    for kk = 1:numel(paths)
        path = paths{kk};
        iy = path(:,1); ix = path(:,2);
        T = ix(end);
        k = numel(ix);
        j = 1;
        for i = 1:T-1
            while(ix(j+1)==i)
                j = j+1;
            end
            midPoints(kk,i) = (iy(j)+iy(j+1))/2;
        end
    end
    midPoints = midPoints -1;
end