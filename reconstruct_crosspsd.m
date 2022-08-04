% read the set
function [S, f] = reconstruct_crosspsd(x, recset, params)
fs = params.fs;
win = params.win;
nov = params.nov;
nfft = params.nfft;
if ~isfield(params, 'method')
    params.method = 'maxdet';
end
compmeth = params.method;

nset = length(recset);

count_mat = meacount_mat(recset);
miss_mat = count_mat == 0;

n = size(count_mat, 1);

for i = 1:nset
    [Si, f] = crosspsd(x{i}, win, nov, nfft, fs);
    nf = length(f);
    if i == 1
        S = zeros(n, n, nf);
    end
    S(recset{i}, recset{i}, :) = S(recset{i}, recset{i}, :) + Si;
end

count_mat(count_mat==0) = 1;
S = S ./ repmat(count_mat, 1, 1, nf);

% determine if completion needed
if ~isempty(find(miss_mat,1))
    disp('matrix incomplete, completing')
    switch compmeth
        case 'maxdet'
            S = maxdet_completion(S, miss_mat);
            disp('max determinant method used')
        case 'mindet'
            S = S;
            disp('min determinant method used')
        otherwise
            disp('Unknown completion method')
    end
else
    % if complete reconstruction, check and correct for positive definit
    S = positivedef_correction(S);
    disp('matrix complete, checking for positive definite')
end