function [F_new,V_new] = MeshSlit(TR,path_list,varargin)
% if 3rd input argument: excludes the start/end nodes of path_list

F = TR.ConnectivityList;
V = TR.Points;
F_new = F;
V_new = V;
E = edges(TR);
NofVertices = size(TR.Points,1);

% Option2: double connected genus 0 mesh
if ~isempty(varargin)
    temp1 = vertexAttachments(TR,path_list);
    temp2 = cellfun(@(x) unique(F(x,:)),temp1,'UniformOutput',0);
    all_neighb = setdiff(unique(cell2mat(temp2(2:end-1))),path_list);
    ind = find(~prod(ismember(E,all_neighb),2));
    Ecropped = E;
    Ecropped(ind,:) = [];
    VconnV = spones(sparse([Ecropped(:,2);Ecropped(:,1)],[Ecropped(:,1);Ecropped(:,2)],[Ecropped(:,2);Ecropped(:,1)],NofVertices,NofVertices));
    [~,C] = graphconncomp(VconnV);
    disc = find(C==min(C(all_neighb)));
    if ~(all(ismember(C(all_neighb),minmax(C(all_neighb)))))
        error('Error in number of connected components');
    end
    
    for k = 2:length(path_list)-1
    adjnodes = intersect(temp2{k},disc);
    current = path_list(k-1:k+1);
    indx = find(sum(ismember(F,[current;adjnodes]),2)==3);
    % add duplicate node
    V_new(end+1,:) = V(current(2),:);
    update_value = size(V_new,1);
    % update 'F' list
    for t = 1:length(indx)
        F_new(indx(t),find(F(indx(t),:)==current(2))) = update_value;
    end
    end
else
    temp1 = vertexAttachments(TR,path_list);
    temp2 = cellfun(@(x) unique(F(x,:)),temp1,'UniformOutput',0);
    all_neighb = setdiff(unique(cell2mat(temp2)),path_list);
    ind = find(~prod(ismember(E,all_neighb),2));
    Ecropped = E;
    Ecropped(ind,:) = [];
    VconnV = spones(sparse([Ecropped(:,2);Ecropped(:,1)],[Ecropped(:,1);Ecropped(:,2)],[Ecropped(:,2);Ecropped(:,1)],NofVertices,NofVertices));
    [~,C] = graphconncomp(VconnV);
    disc = find(C==min(C(all_neighb)));
    disc = disc';
    %****
    if ~(all(ismember(C(all_neighb),minmax(C(all_neighb)))))
        error('Error in number of connected components');
    end
    %****
    for k = 1:length(path_list)
        adjnodes = intersect(temp2{k},[disc;path_list]);
        current = path_list(k);
        indx = find(sum(ismember(F,adjnodes),2)==3);
        % add duplicate node
        V_new(end+1,:) = V(current,:);
        update_value = size(V_new,1);
        % update 'F' list
        for t = 1:length(indx)
            F_new(indx(t),find(F(indx(t),:)==current)) = update_value;
        end
    end
end
%% DIAGNOSTICS
% [boundary1,boundary2] = Detect_Boundaries(TR);
% figure, P = trimesh(TR); axis equal, set(gca,'ZDir','reverse');
% set(P,'FaceColor','interp','FaceVertexCData',cdata,'EdgeColor',[0 0 0],'LineWidth',1,'CDataMapping','scaled')
% for i=1:length(boundary1)
% p1 = plot3(V(boundary1,1),V(boundary1,2),V(boundary1,3),'.r','MarkerSize',20);
% p2 = plot3(V(boundary2,1),V(boundary2,2),V(boundary2,3),'.g','MarkerSize',20);
% end
