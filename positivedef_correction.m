function S = positivedef_correction(S)
[~, n, nf] = size(S);
tol = 1e-5;

lambda = zeros(n, nf);
Us = zeros(n, n, nf);
for i = 1:nf
    [U, T] = schur(S(:,:,i));
    Us(:,:,i) = U;
    lambda(:,i) = diag(T);
end

lambda_nonzero = lambda;
if ~isempty(find(lambda<0,1))
    lambda_nonzero(lambda<0) = 0;
    % smooth
    lambda_sm = smooth(lambda_nonzero(1,:));
    % if still < 0, set to a small number
    lambda_sm(lambda_sm == 0) = tol;
    lambda_nonzero(1,:) = lambda_sm;
    lambda_nonzero = lambda_nonzero ./ repmat(sum(lambda_nonzero), n, 1) .* repmat(sum(lambda), n, 1);
    S_pd = zeros(size(S));
    for i =1:nf
        S1 = Us(:,:,i) * diag(lambda_nonzero(:,i)) * Us(:,:,i)';
        S_pd(:,:,i) = S1;
    end
    S = S_pd;
end