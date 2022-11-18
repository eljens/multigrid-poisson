function A = system_matrix(N)
%SYSTEM_MATRIX
% This function constructs a heptadiagonal sparse banded matrix which
% represents the Laplacian operator in 3D. The matrix is updated to
% impose Neumann conditions on vertical boundaries (East,West,North,
% South) and Dirichlet conditions on the horizontal boundaries (Top and
% Bottom).
%
% Syntax: A = system_matrix(N)
%
% Inputs:
%   N: The dimensions of the coarsest grid. 1 times 3 vector.
%
% Outputs:
%   A: The system matrix representing the discrete Laplacian. 2D sparse matrix.
%
% See also: spdiags
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 18th of November 2022
e = ones(N(1),N(2),N(3));

a = -6*e;

% Improsing the Dirichlet conditions
a(:,:,1) = 1;
a(:,:,N(3)) = 1;
e(:,:,1) = 0;
e(:,:,N(3)) = 0;

bT = e; % Top
bB = e; % Bottom
bW = e; % West
bE = e; % East
bS = e; % South
bN = e; % North

% Imposing the Neumann conditions
% West
bW(1,:,2:end-1) = 0;
bE(1,:,2:end-1) = 2;

% East
bE(N(1),:,2:end-1) = 0;
bW(N(1),:,2:end-1) = 2;

% South
bS(:,1,2:end-1) = 0;
bN(:,1,2:end-1) = 2;

% North
bN(:,N(2),2:end-1) = 0;
bS(:,N(2),2:end-1) = 2;

% Flattening the 3D matrices to 1D
bT = bT(:);
bB = bB(:);
bW = bW(:);
bE = bE(:);
bS = bS(:);
bN = bN(:);
a = a(:);

% Constructing the sparse system matrix with spdiags
A = spdiags([bB([(N(1)*N(2)+1):end 1:N(1)*N(2)]) bS([(N(1)+1):end 1:N(1)]) ...
    bW([2:end 1]) a bE([end 1:(end-1)]) bN([(end-N(1)+1):end 1:end-N(1)]) ...
    bT([(end-N(1)*N(2)+1):end 1:end-N(1)*N(2)])],...
    [-N(1)*N(2),-N(1),-1,0,1,N(1),N(1)*N(2)],prod(N),prod(N));

end