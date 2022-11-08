function A = system_matrix(N)
e = ones(N(1),N(2),N(3));

a = -6*e;

% Dirichlet conditions
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

% Neumann conditions
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

bT = bT(:);
bB = bB(:);
bW = bW(:);
bE = bE(:);
bS = bS(:);
bN = bN(:);
a = a(:);

A = spdiags([bB([(N(1)*N(2)+1):end 1:N(1)*N(2)]) bS([(N(1)+1):end 1:N(1)]) bW([2:end 1]) a...
    bE([end 1:(end-1)]) bN([(end-N(1)+1):end 1:end-N(1)]) bT([(end-N(1)*N(2)+1):end 1:end-N(1)*N(2)])],...
    [-N(1)*N(2),-N(1),-1,0,1,N(1),N(1)*N(2)],prod(N),prod(N));

end