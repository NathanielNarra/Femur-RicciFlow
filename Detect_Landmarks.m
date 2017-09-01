function [FHnode,GTnode,FHboundary,FHsupport] = Detect_Landmarks(TR,cdata)
% LANDMARKS: Detect FH and GT nodes in mesh
%   Detailed explanation goes here

E = edges(TR);
VconnV = sparse([E(:,2);E(:,1)],[E(:,1);E(:,2)],[E(:,2);E(:,1)]);
Umax = max(cdata);
[sorted_values, sorted_index] = sort(cdata,'descend');
initial_thresh = 0.9*Umax;
initial_set = sorted_index(sorted_values>=initial_thresh);
active_node_list = sorted_index(sorted_values<initial_thresh);
adjacent_nodes = setdiff(nonzeros(unique(VconnV(:,initial_set))),initial_set);
current_node = active_node_list(1);

flag = 0;
while ~flag
    if ismember(current_node,adjacent_nodes)
        current_adjacent = setdiff(nonzeros(VconnV(:,current_node)),initial_set);
        if max(cdata(current_adjacent))<=cdata(current_node)
            initial_set = [initial_set;current_node];
        else
            initial_set = [initial_set;current_node];
            flag = 1;
        end
        adjacent_nodes = setdiff(nonzeros(VconnV(:,initial_set)),initial_set);
    end
    active_node_list(1) = [];
    current_node = active_node_list(1);
end

% Calculate FHnode based on surface fitting of the point cloud
x = TR.Points(initial_set,1); 
xc = x-mean(x);
y = TR.Points(initial_set,2); 
yc = y-mean(y);
z = cdata(initial_set);
z = z - min(z);
g = fittype('a*exp(-((x-xm)^2/(2*sx^2) + (y-ym)^2/(2*sy^2)))','dependent','z','independent',{'x','y'},'coefficients',{'a','sx','sy','xm','ym'});
[surffit, ~]= fit([xc(:),yc(:)],z(:),g,'StartPoint',[0.9*max(z),0.1,0.1,0,0]); 
peak = [mean(x) + surffit.xm, mean(y) + surffit.ym];
%figure, plot(surffit,[xc(:),yc(:)],z(:)) %Optional diagnostic
FHnode = nearestNeighbor(TR,peak);
FHsupport = find(pdist2(TR.Points,TR.Points(FHnode,:)) < mean([abs(surffit.sx),abs(surffit.sy)]));
FHboundary = setdiff(nonzeros(unique(VconnV(:,FHsupport))),FHsupport);

if isempty(FHboundary)
    disp('error in calculating FHboundary')
elseif length(FHsupport)<100
    disp('FHsupport length error')
end

sorted_index(ismember(sorted_index,initial_set)) = [];
GTnode = sorted_index(1);
end


