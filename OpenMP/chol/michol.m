function A=michol(A)
% Sobrescriu la part triangular inferior de A amb el factor de Cholesky
%
% Exemple d'us:
%   n=64;
%   A=tril(rand(n));
%   A=A+A'+n*eye(n);
%   b=A*ones(n,1);
%   L=michol(A);
%   L=tril(L);
%   err=norm(ones(n,1)-L'\(L\b))

bs=8;   % block size
[n nn]=size(A);
nb=n/bs;  % number of blocks
for k=1:nb
    A((k-1)*bs+1:k*bs,(k-1)*bs+1:k*bs)=chol(A((k-1)*bs+1:k*bs,(k-1)*bs+1:k*bs))';
    for i=k+1:nb
      A((i-1)*bs+1:i*bs,(k-1)*bs+1:k*bs)=A((i-1)*bs+1:i*bs,(k-1)*bs+1:k*bs)/tril(A((k-1)*bs+1:k*bs,(k-1)*bs+1:k*bs))';
    end
    for i=k+1:nb
        for j=i+1:nb
            A((j-1)*bs+1:j*bs,(i-1)*bs+1:i*bs)=A((j-1)*bs+1:j*bs,(i-1)*bs+1:i*bs)-A((j-1)*bs+1:j*bs,(k-1)*bs+1:k*bs)*A((i-1)*bs+1:i*bs,(k-1)*bs+1:k*bs)';
        end
        A((i-1)*bs+1:i*bs,(i-1)*bs+1:i*bs)=A((i-1)*bs+1:i*bs,(i-1)*bs+1:i*bs)-A((i-1)*bs+1:i*bs,(k-1)*bs+1:k*bs)*A((i-1)*bs+1:i*bs,(k-1)*bs+1:k*bs)';
    end
end
