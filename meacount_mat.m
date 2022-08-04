% count the number of measurements for each pair
% input: sets of simultaneous recording
function count_mat = meacount_mat(recset)
n = max(cellfun(@max, recset));
count_mat = zeros(n);
nset = length(recset);
for i = 1:nset
    inds = nchoosek(recset{i}, 2);
    for j = 1:length(recset{i})
        count_mat(recset{i}(j), recset{i}(j)) = count_mat(recset{i}(j), recset{i}(j)) + 1;
    end
    for j = 1:size(inds,1)
        count_mat(inds(j, 1), inds(j, 2)) = count_mat(inds(j, 1), inds(j, 2)) + 1;
        count_mat(inds(j, 2), inds(j, 1)) = count_mat(inds(j, 2), inds(j, 1)) + 1;
    end
end

