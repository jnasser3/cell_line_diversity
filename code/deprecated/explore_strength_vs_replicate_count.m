function explore_strength_vs_replicate_count(ds)
%Explore the relationship between the number of replicates and the brew
%level signature strength.
%
%Inputs: ds - a structure containing zspc level experimental replicates

TRIALS_PER_REP_COUNT = 10;
num_replicates = size(ds.mat,2);
strength_data = zeros(num_replicates*TRIALS_PER_REPLICATE_COUNT,);

for ii = 1:num_replicates
    for jj = 1:TRIALS_PER_REP_COUNT
        idx = randperm(num_replicates,ii);
        data = ds.mat(:,idx);
        
        strength_data(TRIALS_PER_REP_COUNT*(ii-1) + jj) = [num_replicates,...
            norm(modzs(data),2),...
            norm(mean(data,2),2)];
    end
end



end

