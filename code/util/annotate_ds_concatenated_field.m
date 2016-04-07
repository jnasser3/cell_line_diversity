function ds_out = annotate_ds_concatenated_field(ds,dim,fields,field_name)
%ds_out = annotate_ds_concatenated_field(ds,dim,fields,field_name)
%
%Given a gctstruct, annotates the ds with a concatenation of some set of
%fields already in the struct.
%
%Input:
%       ds: dataset
%       dim: either 'row' or 'column'
%       fields: a cell array of the fields to concatenate
%       field_name: name of the new concatenated field
%
%Output:
%       out_ds: the new annotated data set

%Check that fields exist
if strcmp(dim,'row')
    assert(all(ismember(fields,ds.rhd)),'Some fields not present in ds!');
elseif strcmp(dim,'column')
    assert(all(ismember(fields,ds.chd)),'Some fields not present in ds!');
else
    error('Invalid dimension. Must be either row or column')
end

switch dim
    case 'column'
        idx = ismember(ds.chd,fields);
        desc = ds.cdesc(:,idx);
        rows = mat2cell(desc,ones(size(desc,1),1),size(desc,2));
        rows2 = cellfun(@strjoin,rows,...
            'uni',0);
        
        annot = struct('id',ds.cid,...
            field_name,rows2);
        ds_out = annotate_ds(ds,annot,...
            'dim','column');
        
    case 'row'
        idx = ismember(ds.rhd,fields);
        desc = ds.rdesc(:,idx);
        rows = mat2cell(desc,ones(size(desc,1),1),size(desc,2));
        rows2 = cellfun(@strjoin,rows,...
            'uni',0);
        
        annot = struct('id',ds.rid,...
            field_name,rows2);
        ds_out = annotate_ds(ds,annot,...
            'dim','row');
        
end

end

