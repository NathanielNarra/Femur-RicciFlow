function [TR_flat,TR_3d,v_ix,f_ix] = Path_straighten(TR,TR_3d,FH,GT)
% Modify the mesh so that line between FH-GT contains a sequence of
% edges/nodes. This would act as a shortest path between them.
E = edges(TR);
V3d = TR_3d.Points;
V = TR.Points;
Xp = [V(E(:,1),1),V(E(:,2),1)]';
Yp = [V(E(:,1),2),V(E(:,2),2)]';
Xp(3,:) = NaN;
Yp(3,:) = NaN;
Xp = Xp(:);
Yp = Yp(:);

V_mod = V;
x2 = V([FH,GT],1);
y2 = V([FH,GT],2);

[xi,yi,ii] = polyxpoly(Xp,Yp,x2,y2);
Elist = E';
Elist(3,:) = 0;
Elist = Elist(:);
edge_nodes = [Elist(ii(:,1)),Elist(sum(ii,2))];
for k = 1:size(ii,1)
    temp = edge_nodes(k,:);
    nearest(k,1) = temp(knnsearch(V(temp,:),[xi(k),yi(k)]));
end

%remove start and end nodes locations in xi/yi?
nearest = setdiff(unique(nearest),[FH,GT]);

%find the normal vectors wrt line segment
dx = diff(x2);
dy = diff(y2);
linelength = sqrt(dx.^2+dy.^2);
N1 = [dy/linelength,-dx/linelength];
N2 = [-dy/linelength,dx/linelength];

% construct perpendicular line segments from each nearest node to intersect
% the (FH-GT) line segment
xo = V(nearest,1);
yo = V(nearest,2);
bx = [xo+N1(1),xo+N2(1)]';
bx(3,:) = NaN;
bx = bx(:);
by = [yo+N1(2),yo+N2(2)]';
by(3,:) = NaN;
by = by(:);

% find the new positions of the nearest nodes as intersections between
% perpendicular and line segment
[xii,yii] = polyxpoly(bx,by,x2,y2,'unique');
V_mod(nearest,:) = [xii,yii];

% fix degenerate triangles and self intersections
TR_mod = Mesh_flipedges(triangulation(TR.ConnectivityList,V_mod));
F_mod = TR_mod.ConnectivityList;
Nm = faceNormal(TR_mod);
Nm = Nm(:,3);
if min(Nm)~= 1
    sxnodes = find(Nm==-1);
    for i = 1: length(sxnodes)
        temp = F_mod(sxnodes(i),:);
        mergenodes = temp(ismember(F_mod(sxnodes(i),:),nearest));
        V_mod(mergenodes(2),:)= V_mod(mergenodes(1),:);
    end
end

[tri,bary] = pointLocation(TR,V_mod(nearest,1),V_mod(nearest,2));
V3d(nearest,:) = barycentricToCartesian(TR_3d,tri,bary);

% clean the mesh (remove duplicate nodes)
[Vo,v_ix,f_ix]=unique(V_mod,'rows','stable');

% delete the 2 faces attached to the edge between the two merged nodes
if size(V_mod,1)~=size(Vo,1)
    F_mod(sum(ismember(F_mod,mergenodes),2)==2,:) = [];
end
F = f_ix(F_mod);

TR_flat = triangulation(F,Vo);

% update 3d mesh
V3d = V3d(v_ix,:);
TR_3d = triangulation(F,V3d);

end

