function out = my_pTRPCA_adapter(Y,Xstar,truth,C)
% Adapter for the authors' p-TRPCA implementation.
%
% Xstar and truth are NOT used by the optimization.
% They exist only because all benchmark adapters share one interface.

assert(exist('ptrpca','file') ~= 0, ...
    'Official ptrpca implementation not found.');

assert(isfield(C,'pTRPCA') && C.pTRPCA.finalized, ...
    'p-TRPCA parameters have not been frozen.');

P = C.pTRPCA;

[n1,n2,n3] = size(Y);
nmin = min(n1,n2);

% Balance parameter
lambda = 1/sqrt(max(n1,n2)*n3);

% Paper data-recovery grouping:
% [1:5], [6:10], [11:n_min]
w = zeros(nmin,1);

i1 = min(5,nmin);
i2 = min(10,nmin);

w(1:i1) = P.wValues(1);

if nmin > 5
    w(6:i2) = P.wValues(2);
end

if nmin > 10
    w(11:nmin) = P.wValues(3);
end

opts = struct();
opts.mu       = P.mu;
opts.tol      = P.tol;
opts.rho      = P.rho;
opts.max_iter = P.maxIter;
opts.DEBUG    = 0;

[L,E,obj,err,iter] = ptrpca( ...
    Y,lambda,w,P.p1,P.p2,opts);

assert(isequal(size(L),size(Y)), ...
    'p-TRPCA returned an incompatible tensor size.');

assert(all(isfinite(L(:))), ...
    'p-TRPCA returned nonfinite reconstruction values.');

out = struct();

out.Xhat = L;
out.score = abs(E);

out.extra.obj = obj;
out.extra.err = err;
out.extra.iter = iter;
out.extra.lambda = lambda;
out.extra.wValues = P.wValues;
out.extra.p1 = P.p1;
out.extra.p2 = P.p2;

out.method = 'p-TRPCA';
out.codeSource = 'Yan--Guo authors public implementation';

end