function [disparityMap] = getDisparityMap(midPoints,H,W)
%% alignment result to disparity map
    disparityMap = nan(H,W);
    for k = 1:H
        curMidPoint = midPoints(k,:);
        %%
        x = 0;
        y = 0;
        tempPath = [];
        while(y<W-1)
            ix = [x:floor(curMidPoint(y+1))]';
            iy = ones(numel(ix),1)*y;
            tempPath = [tempPath;ix,iy];
            x = floor(curMidPoint(y+1));
            if(curMidPoint(y+1)~=floor(curMidPoint(y+1)))
                x = x+1;
            end
            y = y+1;
        end
        ix = [x:W-1]';
        iy = ones(numel(ix),1)*(W-1);
        tempPath = [tempPath;ix,iy];
        tempPath = [tempPath(:,2),tempPath(:,1)];
        tempPath = tempPath+1;
        ix = tempPath(:,1);
        iy = tempPath(:,2);

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
end

