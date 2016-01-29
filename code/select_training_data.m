%Takes a ds. Finds all pert/dose/time combinations that appear in at least
%9 cell lines. Then prune to select exactly one of these combinations per
%line. The final product is a set of signatures that appear in each cell
%line exactly once. 

ds = parse_gctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/cp_n5618x978.gctx');
core_lines = parse_grp('/cmap/projects/cell_line_diversity/data/phase1_core_lines.grp');

combo = ds.cdesc(:,ds.cdict('pert_id_dose_time'));

%Find all combos appearing at least 9 times
T = table(unique(combo), countmember(unique(combo), combo),...
    'VariableNames',{'Combo','Count'});
good_idx = T{:,'Count'} >= 9;
good_T = T(good_idx,:);
good_combo = good_T{:,1};

%% For each combo check that it appears in 9 distinct lines and its a cp
great_combo = [];
for ii = 1:numel(good_combo)
    idx = strcmp(ds.cdesc(:,ds.cdict('pert_id_dose_time')),good_combo(ii));
    lines = ds.cdesc(idx,ds.cdict('cell_id'));
    pert_type = unique(ds.cdesc(idx,ds.cdict('pert_type')));
    
    if numel(unique(lines)) >= 9 && strcmp(pert_type,'trt_cp')
        great_combo = [great_combo, good_combo(ii)];
    end
end
great_idx = ismember(ds.cdesc(:,ds.cdict('pert_id_dose_time')),great_combo);

%% For each great combo, only select one of it's signatures in each cell line
for jj = 1:numel(great_combo)
    %get indices of great combo
    idx = strcmp(ds.cdesc(:,ds.cdict('pert_id_dose_time')),great_combo(jj));
    
    med_rep = median(cell2mat(ds.cdesc(idx,ds.cdict('distil_nsample'))));
    
    %find lines in which there is more than one combo
    line_count = countmember(core_lines,ds.cdesc(idx,ds.cdict('cell_id')));
    dup_lines = core_lines(line_count > 1);
    
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

%% Make a gctx
great_ds = ds_slice(ds,'cidx',find(great_idx));
mkgctx('/cmap/projects/cell_line_diversity/data/cmap_training_data/cp_pruned.gctx',great_ds)