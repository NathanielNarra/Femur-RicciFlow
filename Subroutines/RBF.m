function [final_pos, disp] = RBF(V,source,target,scale)
% INPUT
%   V: Nx2 list of node positions/coordinates
%   source: mx2 list of 'm' landmark coordinates in source
%   target: mx2 corresponding landmark positions in target
%   scale: numebr>0 determines the support region 
% OUTPUT
%   disp: Nx2, displacements in each dimension for every node
%   final_pos: Nx2,final positions of nodes in V after RBF based warping

% Calculate coefficients "alpha"
D_lm = squareform(pdist(source));
r = D_lm/scale;
k = ((1-r).^4).*(4*r+1);
k(D_lm>scale) = 0;
alpha = k\(target-source);

% Calculate displacement field
D = pdist2(V,source);
r = D/scale;
K = ((1-r).^4).*(4*r+1);
K(D>scale) = 0;
disp = K*alpha;
final_pos = V+disp;

end

