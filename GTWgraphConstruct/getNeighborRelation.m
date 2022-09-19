function [s,t] = getNeighborRelation(pix,H,W)
% getNeighborRelation: get neighbor relation for each pixel
% a possible edge from s to t 

    [ih0,iw0] = ind2sub([H,W],pix);
    dw = [0,-1];
    dh = [-1,0];
    s = [];
    t = [];
    for k = 1:numel(dw)
        ih = ih0 + dh(k);
        iw = iw0 + dw(k);
        % curValid
        select = ih>0 & iw>0 & ih<=H & iw<=W;
        pixCur = pix(select);
        ih = ih(select);
        iw = iw(select);
        pixN = sub2ind([H,W],ih,iw);
        % in selected region
        select = ismember(pixN,pix);
        s = [s;pixN(select)];
        t = [t;pixCur(select)];
    end
    % map the label
    mapLabel = zeros(H,W);
    mapLabel(pix) = 1:numel(pix);
    s = mapLabel(s);
    t = mapLabel(t);
end