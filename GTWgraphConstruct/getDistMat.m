function d0 = getDistMat(ref,tst,metric)
% getDistMat: Get Euclidean distance between two curves
% ref along height and tst along width

T1 = numel(ref);
T2 = numel(tst);
if(strcmp(metric,'squared'))
    d0 = (repmat(reshape(ref,[],1),1,T2) - repmat(reshape(tst,1,[]),T1,1)).^2;
elseif(strcmp(metric,'absolute'))
    d0 = abs(repmat(reshape(ref,[],1),1,T2) - repmat(reshape(tst,1,[]),T1,1));
end

d0(isnan(d0)) = 0; % !!
d0 = double(d0);

end