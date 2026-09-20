function Xr = hcl_project_tubal_rank(X, r)
% Frobenius projection onto tubal-rank <= r by FFT + per-slice SVD truncation.

[n1,n2,n3] = size(X);
r = min([r,n1,n2]);

Xf = fft(X,[],3);
Zf = zeros(size(Xf),'like',Xf);

for k = 1:n3
    [U,S,V] = svd(Xf(:,:,k),'econ');
    rr = min([r,size(S,1),size(S,2)]);
    Zf(:,:,k) = U(:,1:rr) * S(1:rr,1:rr) * V(:,1:rr)';
end

Xr = real(ifft(Zf,[],3));
end
