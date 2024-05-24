function A = FE(A,tol) % FindEntries
% Converts matrix entries into 0 for entry<tol and 1 for entry>=tol
%INPUT:
%   A = matrix
%   tol = tolerance, below which an entry is considered to be zero
%OUTPUT
%   A = matrix of 0s and 1s such that 0 for entry<tol and 1 for entry>=tol


A = (abs(A)/max(max(abs(A))))>tol;