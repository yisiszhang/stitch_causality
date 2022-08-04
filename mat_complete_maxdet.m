% gd to max determinant
% cost = -det(A)
% grad = -adjoint(A)'
recset = {[1, 2], [1, 3]};
missset = {[2, 3]};

% simple correlation
n = 3;
x0 = eye(n);
recvals = [0.2, 0.3];
for i = 1:length(recset)
    x0(recset{i}(1), recset{i}(2)) = recvals(i);
    x0(recset{i}(2), recset{i}(1)) = recvals(i);
end
ind = ones(size(x0));
% init unknown
for i = 1:length(missset)
    x0(missset{i}(1), missset{i}(2)) = 0.1;
    x0(missset{i}(2), missset{i}(1)) = 0.1;
    ind(missset{i}(1), missset{i}(2)) = 0;
    ind(missset{i}(2), missset{i}(1)) = 0;
end

x0 = [ 0.9874 + 0.0000i  -1.1382 - 0.1101i  -0.2613 - 0.0546i
  -1.1382 + 0.1101i   2.1349 + 0.0000i   0.0000 + 0.0000i
  -0.2613 + 0.0546i   0.0000 + 0.0000i   0.5456 + 0.0000i];

%%
ind = ind == 0;
ind = triu(ind, 1);

%%
niter = 2000;
dt = zeros(1, niter);
alpha = 0.005;
x = x0;

for i = 1:niter
    dt(i) = det(x);
    adjx = adjoint(x)';
    grad = -adjx(ind == 1);
    x(ind == 1) = x(ind == 1) - alpha * grad;
    x(ind' == 1) = conj(x(ind == 1));
end
x1 = x0;
x1(ind == 1) = x0(2,1) * x0(1,3) / x0(1,1);
x1(ind' == 1) = conj(x1(ind == 1));
figure(1)
plot(real(dt))

%%
options = optimset('GradObj', 'on', 'MaxIter', 1000); 
[x, cost] = fminunc(@(x)(costFunction(x, ind)), x0, options);

inv(reshape(x, n, n))

function [cost, grad] = costFunction(x, ind)
    n = length(x);
    x = reshape(x, sqrt(n), sqrt(n));
    cost = - det(x);

    grad = -adjoint(x)';
    grad(ind) = 0;
    grad = grad(:);
end