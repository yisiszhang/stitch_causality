function [S, f] = crosspsd(x, win, nov, nfft, fs)

% x column vectors
[r, c] = size(x);
if c > r
    x = x';
    [r, c] = size(x);
end

% [pxx,f] = cpsd(x, x, win, nov, nfft, fs);
% nf = length(f);
% S = zeros(c, c, nf);
% 
% for i = 1 : c-1
%     S(i, i, :) = pxx(:,i);
%     for j = i+1 : c
%         [pxy,f] = cpsd(x(:,i), x(:,j), win, nov, nfft, fs);
%         S(i, j, :) = pxy;
%         S(j, i, :) = conj(pxy);
%     end
% end
% 
% S(c, c, :) = pxx(:,c);
[pxy,f] = cpsd(x, x, win, nov, nfft,'mimo',fs);
S = permute(pxy, [2,3,1]);