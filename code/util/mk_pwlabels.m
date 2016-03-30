function labels = mk_pwlabels(ids,varargin)
%Given a set of n labels, return an choosek(n,2) set of labels by
%concatenating pairs of labels in the same ordering as produced by tri2vec
%
% mk_pwlabels({'foo','bar','coffee'})
% labels = 
%     'foo_bar'
%     'foo_coffee'
%     'bar_coffee'

params = {'delimiter'};
dflts = {'_'};
args = parse_args(params,dflts,varargin{:});

n = numel(ids);
labels = cell(nchoosek(n,2),1);

kk=1;
for ii = 2:n
    for jj = 1:(ii-1)
        labels(kk) = strcat(ids(jj),args.delimiter,ids(ii));
        kk = kk+1;
    end
end

end
