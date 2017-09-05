function [Y,path_all] = Matching_Slitnodes(DT_template,DT_target,FHnode_template,FHnode_target)

V_template = DT_template.Points;
V_target = DT_target.Points;


%Template
fe = freeBoundary(DT_template);
ConnMatrix = sparse([fe(:,1);fe(:,2)],[fe(:,2);fe(:,1)],1);
GTnode_template = find(ismember(V_template,[0,2*pi],'rows'));
GTnode_template_duplicate = find(ismember(V_template,[0,0],'rows'));
[~, path1, ~] = graphshortestpath(ConnMatrix,GTnode_template,FHnode_template);
[~, path2, ~] = graphshortestpath(ConnMatrix,GTnode_template_duplicate,FHnode_template);
[~, path_all, ~] = graphshortestpath(ConnMatrix,GTnode_template,GTnode_template_duplicate);

% Target
target_bv = freeBoundary(DT_target); % all boundary vertices of target
ConnMatrix_target = sparse([target_bv(:,1);target_bv(:,2)],[target_bv(:,2);target_bv(:,1)],1);
GTnode_target = find(ismember(V_target,[0,max(V_target(V_target(:,1)==0,2))],'rows'));
GTnode_target_duplicate = find(ismember(V_target,[0,0],'rows'));
if isempty(GTnode_target_duplicate) | isempty(GTnode_target)
    GTnode_target = knnsearch(V_target,[0,2*pi]);
    GTnode_target_duplicate = knnsearch(V_target,[0,0]);
end
[~, p1, ~] = graphshortestpath(ConnMatrix_target,GTnode_target,FHnode_target);
[~, p2, ~] = graphshortestpath(ConnMatrix_target,GTnode_target_duplicate,FHnode_target);

Y = V_template(path1,2); 

target_couples = sparse(p1,ones(size(p1)),p2);
slit_bv = unique([p1,p2]);

for i = 2:length(path1(1:end-1))
    
    indx = path2(i);
    point = [0- (1e-6), Y(i)];
    
    [t,bary] = pointLocation(DT_target,point);
    if isnan(t)
        error_path_index = i;
    end
    f = DT_target.ConnectivityList(t,:);
    zero_node = ~(ismember(f,unique(slit_bv)));
    f(zero_node) = [];
    bary(zero_node) = [];
    
    Y(path_all==indx) = sum(V_target(nonzeros(target_couples(f)),2) .* bary');
end
 Y(end+1) = 0;