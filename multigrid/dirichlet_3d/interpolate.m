function uf = interpolate(uc)
%INTERPOLATE
% This function computes a trilinear interpolation of the input matrix.
%
% Syntax: uf = interpolate(uc)
%
% Inputs:
%   uc: The restriction of uf on the coarse domain. 3D Matrix.
%
% Outputs:
%   uf: The interpolation of uc on the fine domain. 3D Matrix.
%
% See also: RESTRICT
%
% Author: Anton Rydahl
% Richard Petersens Plads, bygn. 324 2800 Kgs. Lyngby
% email: rydahlanton@gmail.com
% 28th of October 2022
    N = size(uc);
    uf = zeros(2*N(1)-1,2*N(2)-1,2*N(3)-1);
    uf(1:2:end,1:2:end,1:2:end) = uc;
    uf(2:2:end,:,:) = (uf(3:2:end,:,:) + uf(1:2:end-2,:,:))/2;
    uf(:,2:2:end,:) = (uf(:,3:2:end,:) + uf(:,1:2:end-2,:))/2;
    uf(:,2:2:end,:) = (uf(:,3:2:end,:) + uf(:,1:2:end-2,:))/2;
    uf(:,:,2:2:end) = (uf(:,:,3:2:end) + uf(:,:,1:2:end-2))/2;
end