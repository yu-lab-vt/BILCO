function [midPoints] = path2midPoints(path,T1)
% convert path to midPoint
    ix = path(:,1);
    iy = path(:,2);
    midPoints = zeros(1,T1-1);
    preX = 1;
    k = 2;
    while(k<=numel(ix))
       while(k<=numel(ix) && ix(k)==preX) 
          k = k+1; 
       end
       if(k>numel(ix))
           break;
       end
       midPoints(preX) = (iy(k) + iy(k-1))/2;
       preX = ix(k);
       k = k+1;
    end
end

