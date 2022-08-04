% main script
clear 
close all
% simulate data
% recset = {[1, 2], [1, 3]};
% w = [0; 0; 0];
% A = [0.5 -0.1 0.1; -0.4 0.5 0; -0.1 0.2 0.3];
% C = eye(3)*0.1;

recset = {[1, 2, 3], [2, 3, 4]};
w = [0;0;0;0];
A = [0.5 -0.1 0.1, 0.3; -0.4 0.5 0, 0.2; -0.1 0.2 0.4 0.1; 0.1 0 0.1 0.3];
C = eye(4)*0.1;

n = 10000;
ndisc = 1000;

x = [];
v=arsim(w,A,C,n,ndisc);
x{1} = v(:, recset{1});
v=arsim(w,A,C,n,ndisc);
x{2} = v(:, recset{2});

%%
% reconstruct psd
params.fs = 1;
params.win = bartlett(128);
params.nov = 64;
params.nfft = 1024;
params.method = 'maxdet';
[S, f] = reconstruct_crosspsd(x, recset, params);
disp('Cross psd reconstructed')
%%
c = nonparam_pdc(S, f);
[~, nc, nf] = size(c);
%%
c_truth = zeros(size(c));
for i = 1:nf
 Af = eye(nc) - A * exp(pi * sqrt(-1) / nf * (i-1));
 c_truth(:,:,i) = Af./repmat(sqrt(sum(abs(Af).^2)),nc,1);
end

%%
figure
for i = 1:nc
    for j = 1:nc
        subplot(nc,nc,(i-1)*nc+j)
        plot(squeeze(abs(c(i,j,:))))
        hold on
        plot(squeeze(abs(c_truth(i,j,:))))
    end
end
title('pdc')