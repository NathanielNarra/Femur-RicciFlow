% Similar in size to TR.ConnectivityList. Each row is the Face id and each
% element in the columns contains the Edge id opposite to the vertex in the
% corresposnding position in TR.ConnectivtyList. Edge id's are determined
% by the 'edges' function on TR.

function FE = Faces_constituting_Edges(TR)
F = TR.ConnectivityList;
E = edges(TR);
noE = size(E,1); 
noV = size(TR.Points,1);
% FE = zeros(size(F));

v1 = circshift(F,[0 -1]);
v2 = circshift(F,[0 -2]);
VconnE = full(sparse([E(:,2);E(:,1)],[E(:,1);E(:,2)],[(1:noE)';(1:noE)'],noV,noV));
FE = arrayfun(@(x,y) VconnE(x,y),v1,v2);

% % ***OLD CODE***
% for i = 1:size(F,1)
%     v1 = F(i,1);
%     v2 = F(i,2);
%     v3 = F(i,3);
%     FE(i,1) = find(sum((E==v2 | E==v3),2)==2);
%     FE(i,2) = find(sum((E==v1 | E==v3),2)==2);
%     FE(i,3) = find(sum((E==v1 | E==v2),2)==2);
% end
% % ***