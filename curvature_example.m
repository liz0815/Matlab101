%CURVATURE_EXAMPLE Demonstrate surface curvature computation.
%   This script computes the principal curvatures of the surface
%   z = sin(x) * cos(y) over a rectangular grid and prints the
%   maximum and minimum principal curvature values.

[x,y] = meshgrid(-2:0.1:2);
z = sin(x) .* cos(y);

[k1,k2] = surface_curvature(x,y,z);

maxCurv = max(k1(:));
minCurv = min(k2(:));

fprintf('Max curvature: %.5f\n', maxCurv);
fprintf('Min curvature: %.5f\n', minCurv);
