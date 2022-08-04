function S = maxdet_completion(S, miss_mat)
nf = size(S, 3);

ind = triu(miss_mat, 1);

for i = 1:nf
    S(:,:,i) = gd_maxdet(S(:,:,i), ind);
end
% x0 = S(:);
% 
% options = optimset('GradObj', 'on', 'MaxIter', 1000); 
% [x, cost] = fminunc(@(x)(costFunction(x, ind)), x0, options);
% 
% S = reshape(x, n, n);
% 
% function [cost, grad] = costFunction(x, ind)
%     n = length(x);
%     x = reshape(x, sqrt(n), sqrt(n));
%     cost = - det(x);
% 
%     grad = -adjoint(x)';
%     grad(ind) = 0;
%     grad = grad(:);
% end

end

function x = gd_maxdet(x, ind)

maxiter = 1000;

alpha = 1e-2;
tol = 1e-5;
lastcost = real(-det(x));

for i = 1:maxiter
    % d(det(x))/dx = adjoint(x)'
    adjx = adjoint(x)';
    grad = -adjx(ind == 1);
    x(ind == 1) = x(ind == 1) - alpha * grad;
    x(ind' == 1) = conj(x(ind == 1));
    cost = real(-det(x));
    if lastcost - cost < tol
        break;
    else
        lastcost = cost;
    end
end

end