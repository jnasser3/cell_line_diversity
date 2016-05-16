function out_ds = select_training_data_function(ds,varargin)
%out_ds = select_training_data_function(ds,varargin)
%
%Takes as input an annotated genes x samples struct. Returns a new struct 
%which contains one unique signature per cell line. 'unique' is defined in
%the sense of some column annotation. 
%
%Inputs
%       ds: a genes x samples struct
%       cell_lines: cell array or .grp file. Set of cell id's to consider.
%                   Default: all in ds
%       match_field: field of column annotations to use for matching
%       min_rep: Integer. Minimum number of replicates needed to be
%                      considered as an output signature. Default 2.
%
%Outputs
%       out_ds: a genes x sample struct containing one a unique (in terms
%               of match_field) signature per cell line.
%
%To match uniquely with multiple fields, ie a concatenation of fields,
%combine this function with annotate_ds_concatenated_field.

params = {'cell_lines',...
          'match_field',...
          'min_rep'};
dflts = {'',...
         'pert_id_dose_time',...
         2};
args = parse_args(params,dflts,varargin{:});

cell_lines = get_array_input(args.cell_lines,unique(ds.cdesc(:,ds.cdict('cell_id'))));
num_lines = numel(cell_lines);
match_field = ds.cdesc(:,ds.cdict(args.match_field));

%Subset to cell lines
idx = ismember(ds.cdesc(:,ds.cdict('cell_id')),cell_lines);
ds = ds_slice(ds,'cid',ds.cid(idx));

%remove all signatures not having enough replicates
good_rep_idx = find(cell2mat((ds.cdesc(:,ds.cdict('distil_nsample')))) >= args.min_rep);
ds = ds_slice(ds,'cidx',good_rep_idx);

%Find all combos appearing at least 9 times
T = table(unique(match_field), countmember(unique(match_field), match_field),...
    'VariableNames',{'match_field','Count'});
good_idx = T{:,'Count'} >= num_lines;
good_T = T(good_idx,:);
good_combo = good_T{:,1};

%% For each combo check that it appears in 9 distinct lines
great_combo = good_combo;
for ii = 1:numel(good_combo)
    idx = strcmp(ds.cdesc(:,ds.cdict(args.match_field)),good_combo(ii));
    lines = ds.cdesc(idx,ds.cdict('cell_id'));
    %pert_type = unique(ds.cdesc(idx,ds.cdict('pert_type')));
    
    if numel(unique(lines)) < num_lines %&& strcmp(pert_type,'trt_cp')
        %great_combo = [great_combo, good_combo(ii)];
        great_combo{ii} = 'TOREMOVE';
    end
end
great_combo(strcmp(great_combo,'TOREMOVE')) = [];
great_idx = ismember(ds.cdesc(:,ds.cdict(args.match_field)),great_combo);

%% For each great combo, only select one of it's signatures in each cell line
for jj = 1:numel(great_combo)
    %get indices of great combo
    idx = strcmp(ds.cdesc(:,ds.cdict(args.match_field)),great_combo(jj));
    
    med_rep = median(cell2mat(ds.cdesc(idx,ds.cdict('distil_nsample'))));
    
    %find lines in which there is more than one combo
    line_count = countmember(cell_lines,ds.cdesc(idx,ds.cdict('cell_id')));
    dup_lines = cell_lines(line_count > 1);
    
    %Once we have the duped lines, choose one combo per line
    %Choose the one closest to the median number of replicates
    for kk = 1:numel(dup_lines)
        idx2 = find(idx & strcmp(ds.cdesc(:,ds.cdict('cell_id')),dup_lines(kk)));
         
        num_reps = ds.cdesc(idx2,ds.cdict('distil_nsample'));
        deviation = abs(cell2mat(num_reps) - med_rep);
        [~,min_idx] = min(deviation); %min_idx(1) is the one to keep
        idx_to_keep = idx2(min_idx(1));
        
        %Remove it
        great_idx(idx2) = 0;
        great_idx(idx_to_keep) = 1;  
    end
end

out_ds = ds_slice(ds,'cidx',find(great_idx));

end