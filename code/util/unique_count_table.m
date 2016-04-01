function T = unique_count_table(ds,dim,field)
%T = unique_count_table(ds,dim,field)
%
%Make a table summarizing the number of times an annotation appears
%in a ds
%
%Input:
%       ds - An annotated struct
%       dim - Either 'row' or 'column'
%       field - the field name of the annotation
%
%Output:
%       T - A table whose row names contain the unique annotations. T has
%       one column giving the number of times each annotation appears in
%       the dataset

switch dim
    case 'row'
        annots = ds.rdesc(:,ds.rdict(field));
    case 'column'
        annots = ds.cdesc(:,ds.cdict(field));
end

counts = countmember(unique(annots),annots);
T = array2table(counts,'RowNames',unique(annots));


end

