function tDly = getDelayMidPoints(midPoints,ref,thrVec)
% getDelay: estimate the rising time
% INPUT
% paths: warping functions
% ref:   reference curve
% thrVec: tested threshold hold

% OUTPUT
% tDly: estimated rising time for each curve

    [nSp,T] = size(ref);
    paths = cell(nSp,1);
    
    for i = 1:nSp
        curMidPoint = midPoints(i,:);
        %%
        x = 0;
        y = 0;
        tempPath = [];
        while(y<T-1)
            ix = [x:floor(curMidPoint(y+1))]';
            iy = ones(numel(ix),1)*y;
            tempPath = [tempPath;ix,iy];
            x = floor(curMidPoint(y+1));
            if(curMidPoint(y+1)~=floor(curMidPoint(y+1)))
                x = x+1;
            end
            y = y+1;
        end
        ix = [x:T-1]';
        iy = ones(numel(ix),1)*(T-1);
        tempPath = [tempPath;ix,iy];
        tempPath = [tempPath(:,2),tempPath(:,1)];
        paths{i} = tempPath+1;
    end
    
    cx = zeros(nSp,T);
    for i = 1:nSp
        cx(i,:) = warpRef2Tst(paths{i},ref(i,:));
    end

    tAch = nan(nSp,numel(thrVec));
    for nn=1:nSp
        x = cx(nn,:);
        [maxV,t0] = max(x);
        x = x(1:t0);
        for ii=1:numel(thrVec)
            t1 = find(x>=maxV*thrVec(ii),1);
            if isempty(t1)
                t1 = t0;
            end
            tAch(nn,ii) = t1;
        end
    end
    tDly = mean(tAch,2);
end