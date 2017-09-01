% x,y,z are Nx1 vectors. coordinates of N vertices
% TRI is a Mx3 array where M is number of triangles/faces
% Output matrix: [ELjk ELik ELij] where EL is edge length and i,j,k denote
% vertices of the triangular face. Mx3 where M is the number of faces
function [edge_length, edge_list_length] = Discrete_Riemannian_Metric(TR)
%% OLD version (26 seconds for 50k mesh)
% for i = 1:size(F,1)
%     for j = 1:3
%         mesh(i,j,:) = [x(F(i,j)) y(F(i,j)) z(F(i,j))];
%     end
% end
% 
% mesh_shift = circshift(mesh,[0 -1 0]);
% edge_length = sqrt(sum((mesh-mesh_shift).^2,3));
% edge_length = circshift(edge_length,[0 -1]);
% 
% % Edge lengths
% E = edges(TR);
% t1 = num2cell(TR.Points,2);
% %t2 = t1(E);
% t = permute(cell2mat(permute(t1(E),[1 3 2])),[1 3 2]);
% edge_list_length = sqrt(sum((t(:,1,:) - t(:,2,:)).^2,3));

%% OPTIMIZED 1 (10 seconds for 50k mesh)
% A = arrayfun(@(x) {TR.Points(x,:)},TR.ConnectivityList);
% v1 = circshift(A,[0 -1]);
% v2 = circshift(A,[0 -2]);
% edge_length = cell2mat(cellfun(@(x,y) sqrt(sum((x-y).^2)),v1,v2,'UniformOutput',0));

%% OPTIMIZED 2 (<2 seconds for 50k mesh)
E = edges(TR);
V = TR.Points;
v1 = V(E(:,1),:);
v2 = V(E(:,2),:);
edge_list_length = sqrt(sum((v1-v2).^2,2));
FE = Faces_constituting_Edges(TR);
edge_length = edge_list_length(FE);
