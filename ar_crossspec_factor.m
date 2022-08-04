% simulate a 3 way AR model
clear
close all
rng(1)
w = [0; 0; 0];
A = [0.5 -0.1 0.1; -0.4 0.5 0; -0.1 0.2 0.3];
% A = eye(3) * 0.3 + randn(3)*0.01;
%A = eye(3) * 0.3;
%A(1,2) = -0.2;
C = eye(3)*0.1;%[0.1 0.01 0.3; 0.01 0.5 0; 0.3 0 1];

%w = [0;0;0;0];
%A = [0.5 -0.1 0.1, 0.3; -0.4 0.5 0, 0.2; -0.1 0.2 0.4 0.1; 0.1 0 0.1 0.3];
%C = eye(4)*0.1;
n = 10000;
ndisc = 1000;
v=arsim(w,A,C,n,ndisc);

%%
% reconstruct the cross covariance? or cross power spectrum
fs = 1;
win = bartlett(128);
nov = 64;
nfft = 2048;
[pxy,f] = cpsd(v, v, win, nov, nfft,'mimo',fs);

% % %%
% % cross psd using chronux
% win = 1000; 
% params.tapers = [6.5, 12]; 
% [pxy,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(v,win,params);
%%
% spectral factorization
S = permute(pxy, [2,3,1]);
[H, Z, ~, psi] = sfactorization_wilson(S,f);

%%
% simulate pairwise
recsets = {[1,2], [2,3], [1,3]};
pxy_construct = zeros(size(pxy));
for i = 1:length(recsets)
    vr=arsim(w,A,C,n,ndisc);
%     vr =v;
    [pxy1,f] = cpsd(vr(:, recsets{i}(1)), vr(:, recsets{i}(2)), win, nov, nfft, fs);
    pxy_construct(:, recsets{i}(1), recsets{i}(2)) = pxy1;
    pxy_construct(:, recsets{i}(2), recsets{i}(1)) = conj(pxy1);
    
    [pxx,f] = cpsd(vr(:, recsets{i}(1)), vr(:, recsets{i}(1)), win, nov, nfft, fs);
    pxy_construct(:, recsets{i}(1), recsets{i}(1)) = pxy_construct(:, recsets{i}(1), recsets{i}(1)) + abs(pxx)/2;
    [pxx,f] = cpsd(vr(:, recsets{i}(2)), vr(:, recsets{i}(2)), win, nov, nfft, fs);
    pxy_construct(:, recsets{i}(2), recsets{i}(2)) = pxy_construct(:, recsets{i}(2), recsets{i}(2)) + abs(pxx)/2;
end

% % %%
% % simulate pairwise using chronux
% recsets = {[1,2], [2,3], [1,3]};
% pxy_construct = zeros(size(pxy));
% for i = 1:length(recsets)
%     vr=arsim(w,A,C,n,ndisc);
% %     vr =v;
%     [pxy1,Cmat,Ctot,Cvec,Cent,f]=CrossSpecMatc(vr(:, recsets{i}),win,params);
%     pxy_construct(:, recsets{i}(1), recsets{i}(2)) = pxy1(:, 1, 2);
%     pxy_construct(:, recsets{i}(2), recsets{i}(1)) = conj(pxy1(:, 1, 2));
% 
%     pxy_construct(:, recsets{i}(1), recsets{i}(1)) = pxy_construct(:, recsets{i}(1), recsets{i}(1)) + abs(pxy1(:, 1, 1))/2;
%     pxy_construct(:, recsets{i}(2), recsets{i}(2)) = pxy_construct(:, recsets{i}(2), recsets{i}(2)) + abs(pxy1(:, 2, 2))/2;
% end
% %%
% % simulate pairwise
% recsets = {[1,2], [2,3], [1,3]};
% pxy_construct = zeros(size(pxy));
% for i = 1:length(recsets)
%     vr=arsim(w,A,C,n,ndisc);
%     
%     [pxy1,f] = cpsd(vr(:, recsets{1}(1)), vr(:, recsets{1}(2)), win, nov, nfft, fs);
%     pxy_construct(:, recsets{i}(1), recsets{i}(2)) = pxy1;
%     pxy_construct(:, recsets{i}(2), recsets{i}(1)) = conj(pxy1);
%     
%     [pxx,f] = cpsd(vr(:, recsets{i}(1)), vr(:, recsets{i}(1)), win, nov, nfft, fs);
%     pxy_construct(:, recsets{i}(1), recsets{i}(1)) = pxx;
%     [pxx,f] = cpsd(vr(:, recsets{i}(2)), vr(:, recsets{i}(2)), win, nov, nfft, fs);
%     pxy_construct(:, recsets{i}(2), recsets{i}(2)) = pxx;
% 
%         
% end

%%
S_construct = permute(pxy_construct, [2,3,1]);
% [H_construct, Z_construct, ~, psi_construct] = sfactorization_wilson(S_construct,f);

%%
nf = size(S_construct,3);
lambda = zeros(3, nf);
Us = zeros(3, 3, nf);
for i = 1:nf
    [U, T] = schur(S_construct(:,:,i));
    Us(:,:,i) = U;
    lambda(:,i) = diag(T);
end

%%
lambda_nonzero = lambda;
if ~isempty(find(lambda<0))
    lambda_nonzero(lambda<0) = 0;
    % smooth till no 0
    lambda_sm = smooth(lambda_nonzero(1,:));
    lambda_sm(lambda_sm == 0) = 1e-5;
    lambda_nonzero(1,:) = lambda_sm;
    lambda_nonzero = lambda_nonzero ./ repmat(sum(lambda_nonzero), 3, 1) .* repmat(sum(lambda), 3, 1);
else
    lambda_nonzero = lambda;
end
S_construct_pd = zeros(size(S_construct));

for i =1:nf
    S1 = Us(:,:,i) * diag(lambda_nonzero(:,i)) * Us(:,:,i)';
    S_construct_pd(:,:,i) = S1;
end
%%
[H_construct, Z_construct, ~, psi_construct] = sfactorization_wilson(S_construct_pd,f);
%%
Hf = zeros(3,3,nf);
for i = 1:nf
 Hf(:,:,i) = inv(eye(3) - A * exp(pi * sqrt(-1) / nf * (i-1)));
end

%%
figure
plot(squeeze(abs(H(2,3,:))))
hold on
plot(squeeze(abs(H_construct(2,3,:))))
plot(squeeze(abs(Hf(2,3,:))))
title('H')
%%
lambda_simul = zeros(3, nf);
for i = 1:nf
    lambda_simul(:,i) = eig(S(:,:,i));
end

figure
subplot(121)
plot(lambda_simul')
title('eigenvalue simultaneous')
subplot(122)
plot(lambda')
title('eigenvalue reconstruct')
%%
% pdc
pdc = zeros(size(H_construct));
for i = 1:nf
    Af = inv(H_construct(:,:,i));
    pdc(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),3,1);
end

%%
pdc_real = zeros(3,3,nf);
for i = 1:nf
 Af = eye(3) - A * exp(pi * sqrt(-1) / nf * (i-1));
 pdc_real(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),3,1);
end

%%
pdc_simul = zeros(3,3,nf);
for i = 1:nf
 Af = inv(H(:,:,i));
 pdc_simul(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),3,1);
end
%%
figure
for i = 1:3
    for j = 1:3
        subplot(3,3,(i-1)*3+j)
        plot(squeeze(abs(pdc(i,j,:))))
        hold on
        plot(squeeze(abs(pdc_simul(i,j,:))))
        plot(squeeze(abs(pdc_real(i,j,:))))
    end
end
title('pdc')
%%

% simulate missing one pair with min determinant
recsets = {[1,2], [1,3]};
pxy_construct_incom = NaN(size(pxy));
for i = 1:length(recsets)
    vr=arsim(w,A,C,n,ndisc);
%     vr =v;
    [pxy1,f] = cpsd(vr(:, recsets{i}(1)), vr(:, recsets{i}(2)), win, nov, nfft, fs);
    pxy_construct_incom(:, recsets{i}(1), recsets{i}(2)) = pxy1;
    pxy_construct_incom(:, recsets{i}(2), recsets{i}(1)) = conj(pxy1);
    
    [pxx,f] = cpsd(vr(:, recsets{i}(1)), vr(:, recsets{i}(1)), win, nov, nfft, fs);
    pxy_construct_incom(:, recsets{i}(1), recsets{i}(1)) = abs(pxx);
    [pxx,f] = cpsd(vr(:, recsets{i}(2)), vr(:, recsets{i}(2)), win, nov, nfft, fs);
    pxy_construct_incom(:, recsets{i}(2), recsets{i}(2)) = abs(pxx);
end

%%
nf = length(f);
pxy_construct_incom_maxdet = pxy_construct_incom;
for i = 1:nf
    p1 = squeeze(pxy_construct_incom(i,:,:));
    pxy_construct_incom(i,2,3) = 0;
    pxy_construct_incom(i,3,2) = 0;

    x = p1(2,1) * p1(1,3) / p1(1,1);
    pxy_construct_incom_maxdet(i,2,3) = x;
    pxy_construct_incom_maxdet(i,3,2) = conj(x);
end

%%
S_construct_incom = permute(pxy_construct_incom, [2,3,1]);
S_construct_incom_maxdet = permute(pxy_construct_incom_maxdet, [2,3,1]);
%%
nf = size(S_construct_incom,3);
lambda = zeros(3, nf);

for i = 1:nf
    [U, T] = schur(S_construct_incom(:,:,i));
    lambda(:,i) = diag(T);
end

%%
lambda_maxdet = zeros(3, nf);

for i = 1:nf
    [U, T] = schur(S_construct_incom(:,:,i));
    lambda_maxdet(:,i) = diag(T);
end
%%
[H_construct, Z_construct, ~, psi_construct] = sfactorization_wilson(S_construct_incom,f);
[H_construct_maxdet, Z_construct, ~, psi_construct] = sfactorization_wilson(S_construct_incom_maxdet,f);
% %%
% Hf = zeros(3,3,nf);
% for i = 1:nf
%  Hf(:,:,i) = inv(eye(3) - A * exp(pi * sqrt(-1) / nf * (i-1)));
% end
% 
% %%
% figure
% plot(squeeze(abs(H(2,3,:))))
% hold on
% plot(squeeze(abs(H_construct(2,3,:))))
% plot(squeeze(abs(Hf(2,3,:))))
% title('H')
%%
lambda_simul = zeros(3, nf);
for i = 1:nf
    lambda_simul(:,i) = eig(S(:,:,i));
end

figure
subplot(131)
plot(lambda_simul')
title('eigenvalue simultaneous')
subplot(132)
plot(lambda')
title('eigenvalue reconstruct mindet')
subplot(133)
plot(lambda_maxdet')
title('eigenvalue reconstruct maxdet')
%%
% pdc
pdc_mindet = zeros(size(H_construct));
for i = 1:nf
    Af = inv(H_construct(:,:,i));
    pdc_mindet(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),3,1);
end

%%
% pdc
pdc_maxdet = zeros(size(H_construct_maxdet));
for i = 1:nf
    Af = inv(H_construct_maxdet(:,:,i));
    pdc_maxdet(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),3,1);
end

%%
pdc_real = zeros(3,3,nf);
for i = 1:nf
 Af = eye(3) - A * exp(pi * sqrt(-1) / nf * (i-1));
 pdc_real(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),3,1);
end

%%
pdc_simul = zeros(3,3,nf);
for i = 1:nf
 Af = inv(H(:,:,i));
 pdc_simul(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),3,1);
end
%%
figure
for i = 1:3
    for j = 1:3
        subplot(3,3,(i-1)*3+j)
        plot(squeeze(abs(pdc_mindet(i,j,:))))
        hold on
        plot(squeeze(abs(pdc_maxdet(i,j,:))))
        plot(squeeze(abs(pdc(i,j,:))))
        plot(squeeze(abs(pdc_real(i,j,:))))
    end
end
title('pdc')
