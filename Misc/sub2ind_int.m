function ndx = sub2ind_int(siz,v1,v2,v3)
%   SUB2IND_INT computes the linear index from subvectors
%
%   ndx = sub2ind_int(siz, v1, v2, v3)
%   
%   This function works the same as sub2ind, but works for integers,
%   resulting in a faster construction.

    k = int32(cumprod(siz));
    ndx = v1 + k(1)*(v2-int32(1)) + k(2)*(v3-int32(1));
end