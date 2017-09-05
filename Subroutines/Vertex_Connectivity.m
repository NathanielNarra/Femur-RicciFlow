% Connectivity Matrix
% INPUT: {TR} Triangulation of a mesh
% OUTPUT: 
% V_conn_matrix is sparse matrix 
function [V_conn_matrix, V_conn_matrix_mask, V_ind_list, V_ind_mask]= Vertex_Connectivity(TR)

F = TR.ConnectivityList;
E = edges(TR);
noe = size(E,1);
nov = size(TR.Points,1);
nof = size(TR.ConnectivityList,1);
V_conn_matrix = sparse([E(:,1);E(:,2)],[E(:,2);E(:,1)],[(1:noe)';(1:noe)'],nov,nov);
index_image = bsxfun(@plus,repmat((1:nof)',[1 3]),[0,nof,2*nof]);
V_ind_list = sparse(repmat((1:nof)',[3 1]),F(:),index_image(:));
V_conn_matrix_mask = spones(V_conn_matrix);
V_ind_mask = spones(V_ind_list);
