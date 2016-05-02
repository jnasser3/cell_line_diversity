function ds = ds_subset(ds,dim,field,to_include)
% ds = ds_subset(ds,dim,field,to_include)
% 
% Subset a ds based on meta data.

    switch dim
        case 'column'
            good_idx = find(ismember(ds.cdesc(:,ds.cdict(field)),to_include));
            ds = ds_slice(ds, 'cidx', good_idx);
        case 'row'
            good_idx = find(ismember(ds.rdesc(:,ds.rdict(field)),to_include));
            ds = ds_slice(ds, 'ridx', good_idx);
    end
end
