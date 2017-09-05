function [ TR_mod ] = Mesh_flipedges( TR )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

F = TR.ConnectivityList;
V = TR.Points;
[EL, ~] = Discrete_Riemannian_Metric(TR);
%FE = Faces_constituting_Edges(TR);
numerator = repmat(sum(EL.^2,2),[1,3]) - 2*(EL.^2);
denominator = repmat(2*prod(EL,2),[1,3])./EL;
Q = numerator./denominator;
if min(Q(:)) < -0.98
    
    ix = Q<-0.98;
    %Vix = F(ix);
    Tix = find(any(ix,2));
    %flags = ix(Tix,:);
    for i = 1:length(Tix)
        flags = ix(Tix(i),:);
        node_seq = F(Tix(i,:),:);
        edge_nodes = node_seq(~flags);
        % find triangles adjacent to edge
        Tn = cell2mat(edgeAttachments(TR,edge_nodes));
        all4nodes = unique(F(Tn,:));
        %really laborious method
        T1 = F(Tn(1),:);
        toinsert = setdiff(all4nodes,T1);
        toreplace = edge_nodes(1);
        T1(T1==toreplace) = toinsert;
        T2 = F(Tn(2),:);
        toinsert = setdiff(all4nodes,T2);
        toreplace = edge_nodes(2);
        T2(T2==toreplace) = toinsert;
        F(Tn(1),:) = T1;
        F(Tn(2),:) = T2;
        
        %         % images
        %         figure, hold on
        %         plot3(V(edge_nodes,1),V(edge_nodes,2),V(edge_nodes,3))
        %         plot3(V(edge_nodes,1),V(edge_nodes,2),V(edge_nodes,3))
        %         plot3(V(T1(1),1),V(T1(1),2),V(T1(1),3),'.r','MarkerSize',16)
        %         plot3(V(T1(2),1),V(T1(2),2),V(T1(2),3),'.g','MarkerSize',16)
        %         plot3(V(T1(3),1),V(T1(3),2),V(T1(3),3),'.b','MarkerSize',16)
        %         plot3(V(T2(1),1),V(T2(1),2),V(T2(1),3),'or','MarkerSize',16)
        %         plot3(V(T2(2),1),V(T2(2),2),V(T2(2),3),'og','MarkerSize',16)
        %         plot3(V(T2(3),1),V(T2(3),2),V(T2(3),3),'ob','MarkerSize',16)
        %         axis equal
        
    end
end
TR_mod = triangulation(F,V);
end

