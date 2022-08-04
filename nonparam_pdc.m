function c = nonparam_pdc(S, f)
[~, n, nf] = size(S);

% Wilson decomposition
[H, ~, ~, ~] = sfactorization_wilson(S,f);

c = zeros(size(H));
for i = 1:nf
    Af = inv(H(:,:,i));
    c(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),n,1);
end