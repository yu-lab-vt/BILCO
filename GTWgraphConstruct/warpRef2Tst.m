function [x0] = warpRef2Tst(path,ref)
%   warpRef2Tst: warp the reference curve to each pixel
%   Some other features:
%   Generate the new position of given list of time points
%   Input path in 2D cell array, but allow two output types
%   Allow missing values in reference curve

    T = numel(ref);
    x0 = nan(1,T);  % warped curve
    c0 = zeros(1,T);  % count the occurrence
    p0 = path;
    idxValid = p0(:,1)>=1 & p0(:,1)<=T & p0(:,2)>=1 & p0(:,2)<=T;
    p0 = p0(idxValid,:);
    for tt=1:length(p0)
        p_ref = p0(tt,1);
        p_tst = p0(tt,2);
        if ~isnan(ref(p_ref))
            if isnan(x0(p_tst))
                x0(p_tst) = ref(p_ref);
            else
                x0(p_tst) = x0(p_tst) + ref(p_ref);
            end
            c0(p_tst) = c0(p_tst) + 1;
        end
    end
    
    c0(c0==0) = 1;
    x0 = x0./c0;
end