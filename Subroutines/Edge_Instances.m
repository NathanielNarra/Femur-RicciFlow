% returns a matrix listing the indices of the Edges that are in FE (matrix
% of 3 edges that make each face)
function EF = Edge_Instances(TR,FE)

nof = size(TR.ConnectivityList,1);
noe = max(FE(:));
index_image = bsxfun(@plus,repmat((1:nof)',[1 3]),[0,nof,2*nof]);
EF = sparse(repmat((1:nof)',[3 1]),FE(:),index_image(:),nof,noe);
