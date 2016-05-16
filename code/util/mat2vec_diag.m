function [diag_elements, off_diag] = mat2vec_diag(X)
%[diag_elements, off_diag] = mat2vec_diag(X)
%
%Decomposes a square matrix into it's diagonal and off diagonal elements
%
%Note that if X is symmetric, the upper and lower triangles of X are the
%same and thus off_diag is redundant. In this case it may be better to use
%off_diag = tri2vec(X).
%
%Input:
%       X: A square matrix
%
%Output: 
%       diag_elements: A vector containing the diagonal elements of X
%       off_diag: A vector containing the off diagonal elements of X.

[n,m] = size(X);
assert(n == m, 'Input matrix must be square')

diag_elements = diag(X);
off_diag = [tri2vec(X,1,true); tri2vec(X,1,false)];


end

