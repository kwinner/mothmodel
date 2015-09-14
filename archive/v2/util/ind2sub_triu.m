function [ x, y ] = ind2sub_triu( M, i )
%IND2SUB_TRIU similar to ind2sub but for linear indices into the upper triangular part of a square matrix
%   M - side length of the (square) matrix
%   i - the index(-ices) into an upper tri. matrix of size M x M

%N is the number of elements in the upper triangular part of an M x M matrix
N = M*(M+1)/2;

%k is the largest size of a sq.tri. matrix which could hold all elements in rows after the one of i
k = floor(.5*(-1+sqrt(1+8*(N-i))));

%x is the row before the kxk submatrix
x = M-k;

%get the offset from the start of the row by subtracting the starting position of the row of i, then correct for the start of the row
y = i - N + (k+1).*(k+2)/2 + x - 1;

end