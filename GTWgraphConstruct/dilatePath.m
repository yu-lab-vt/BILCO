function [lB,uB] = dilatePath(p,maxShift)
% dilatePath: allow shift based on the esitmated path.
% INPUT
% p: estimated path
% maxShift: allowed gamma

% OUTPUT
% lB: lower bound path
% uB: upper bound path

    if(isempty(p))
        lB = [];
        uB = [];
        return;
    end
    
    T1 = p(end,1);
    T2 = p(end,2);
    p = path2midPoints(p,T1)-1;
    if(maxShift<=0)
        lB = p;
        uB = p;
        return;
    end
    if(maxShift>=max(T1,T2))
        lB = zeros(1,T1-1);
        uB = (T2-1)*ones(1,T1-1);
        return;
    end
    %% lower
    lB = max(0,p-maxShift);
    lB(maxShift+1:end) = min(lB(maxShift+1:end),p(1:end-maxShift));
    lB(1:maxShift) = 0;
    %% upper
    uB = min(T2-1,p+maxShift);
    uB(1:end-maxShift) = max(p(maxShift+1:end),uB(1:end-maxShift));
    uB(end+1-maxShift:end) = T2-1;
    %% 
    jump = floor(p(2:end)) - ceil(p(1:end-1));
    jumpPoint = find(jump>0)+1;
    for i = 1:numel(jumpPoint)
        curP = jumpPoint(i);
        t1 = min(T1 - 1,curP + maxShift-1);
        lB(curP:t1) = min(lB(curP:t1),max(0,ceil(p(curP-1)) - maxShift + 0.5 + [0:t1-curP]));
        t0 = max(1,curP-maxShift);
        uB(t0:curP-1) = max(uB(t0:curP-1),min(T2-1,floor(p(curP)) + maxShift - 0.5 + [t0-curP+1:0]));
    end
end