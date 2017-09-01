function [ ring ] = Ring_Neighbours(TR,N,varargin)
%List of neighbouring nodes to each node in 1...N rings
% 'vargin': number or vector to choose which node(s) to find neighbours for
NofVertices = size(TR.Points,1);
E = edges(TR);
VconnV = sparse([E(:,1);E(:,2)],[E(:,2);E(:,1)],[E(:,1);E(:,2)]);

if nargin==2
    ring = num2cell(1:NofVertices);
    for i = 1:N
        ring(i+1,:) = cellfun(@(x) setdiff(nonzeros(unique(VconnV(:,x))),x),ring(i,:),'UniformOutput',false);
    end
elseif nargin==3
    ring = num2cell(varargin{1});
    for i = 1:N
        ring(i+1,:) = cellfun(@(x) setdiff(nonzeros(unique(VconnV(:,x))),x),ring(i,:),'UniformOutput',false);
    end
end


end

