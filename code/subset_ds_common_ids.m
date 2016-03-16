function [ds1, ds2, common_ids] = subset_ds_common_ids(ds1, ds2, varargin )
%
% Subset to ds's to common ids

params = {'ds1_id_field',...
          'ds2_id_field',...
          'ids'};
dflts = {'cid',...
         'cid',...
         ''};
args = parse_args(params,dflts,varargin{:});

% Get the ds from gctx if necessary
if ischar(ds1)    
    ds1 = parse_gctx(ds1);
end
if ischar(ds2)   
    ds2 = parse_gctx(ds2);
end

%make the cid in each ds equal to the common identifier
if ~strcmp(args.ds1_id_field,'cid')
    ds1.cid = ds1.cdesc(:,ds1.cdict(args.ds1_id_field));
    ds1.rid = ds1.cid;
end
if ~strcmp(args.ds2_id_field,'cid')
    ds2.cid = ds2.cdesc(:,ds2.cdict(args.ds2_id_field));
    ds2.rid = ds2.cid;
end

%get the cell lines to include
common_ids = intersect(ds1.cid, ds2.cid);
if isempty(args.ids)
    args.ids = common_ids;
else
    try
        args.ids = parse_grp(args.ids);
    end
end

%Subset each ds to common lines
%Note that ds_slice will order the structs alphabetically.
%common_lines = intersect(feature_ds.cid, diversity_ds.cid);
common_ids = intersect(common_ids,args.ids);
ds1 = ds_slice(ds1,'rid',common_ids,'cid',common_ids);
ds2 = ds_slice(ds2,'rid',common_ids,'cid',common_ids);

fprintf('ds1 contains %d objects \n',numel(ds1.cid));
fprintf('ds2 contains %d objects \n',numel(ds2.cid));
fprintf('There are %d objects in common \n',numel(common_ids));


end

