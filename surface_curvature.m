function [k1, k2, H, K] = surface_curvature(X, Y, Z)
%SURFACE_CURVATURE Compute principal curvatures of a surface.
%   [K1, K2, H, K] = SURFACE_CURVATURE(X,Y,Z) computes the principal
%   curvatures K1 and K2, mean curvature H, and Gaussian curvature K of the
%   surface defined by arrays X, Y and Z (as returned by MESHGRID). The
%   method uses finite differences to approximate derivatives.
%
%   Example:
%       [x,y] = meshgrid(-1:0.05:1);
%       z = x.^2 + y.^2; % paraboloid
%       [k1,k2,H,K] = surface_curvature(x,y,z);
%       maxCurv = max(k1(:));
%       minCurv = min(k2(:));
%
%   References:
%       do Carmo, "Differential Geometry of Curves and Surfaces".
%
%   See also GRADIENT.

% Determine grid spacing (assumes uniform grid)
if numel(unique(diff(X(1,:)))) > 1 || numel(unique(diff(Y(:,1)))) > 1
    error('X and Y must be uniform grids as produced by MESHGRID.');
end

dx = X(1,2) - X(1,1);
% Y varies along rows
dy = Y(2,1) - Y(1,1);

% First derivatives
[fx, fy] = gradient(Z, dx, dy);

% Second derivatives
[fxx, fxy] = gradient(fx, dx, dy);
[~, fyy] = gradient(fy, dx, dy);

% First fundamental form coefficients
E = 1 + fx.^2;
F = fx .* fy;
G = 1 + fy.^2;

% Normalization factor
W = sqrt(1 + fx.^2 + fy.^2);

% Second fundamental form coefficients
L = fxx ./ W;
M = fxy ./ W;
N = fyy ./ W;

% Mean and Gaussian curvature
H = (E .* N - 2 * F .* M + G .* L) ./ (2 * (E .* G - F.^2));
K = (L .* N - M.^2) ./ (E .* G - F.^2);

% Principal curvatures
Delta = sqrt(max(H.^2 - K, 0));
k1 = H + Delta;
k2 = H - Delta;
end
