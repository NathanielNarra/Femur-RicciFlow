function Flat_V = Seed_face(EL,TR)

Faces = TR.ConnectivityList;
seed = 1000; % embedding starts with triangle face# 1000 (random choice)

% Euclidean geometry: embedding the first face onto a plane
Li = EL(seed,1); % Ljk
Lj = EL(seed,2); % Lik
Lk = EL(seed,3); % Lij
theta =  acos((Lj^2 + Lk^2 - Li^2)/(2*Lj*Lk));
V1 = [0 0];
V2 = [Lk 0];
V3 = Lj*[cos(theta) sin(theta)];

% Assigning vertex indices flags to indicate embedding
Flags = Faces;
Flags(Flags==Faces(seed,1)) = 0;
Flags(Flags==Faces(seed,2)) = 0;
Flags(Flags==Faces(seed,3)) = 0;

%Flat_V = TR.Points;
Flat_V(Faces(seed,1),:) = V1;
Flat_V(Faces(seed,2),:) = V2;  
Flat_V(Faces(seed,3),:) = V3;

% Seed face is the first face from 'TR.ConnectivityList'
Face_id = seed;
Que = [];
Face_id_counter = zeros(1,size(Faces,1));
Face_id_counter(Face_id) = 1;

while nnz(Flags)>0 %do while Queue is populated
    
    %Find neighbouring faces wrt to current face_id.
    nbrs = neighbors(TR,Face_id);
    nbrs = nbrs(~isnan(nbrs)); % removing NaN's
  
    % Append neighbours to end of list
    Que = [Que nbrs];
    Face_id = Que(1);
    Que(1) = [];
    
    while nnz(Flags(Face_id,:))==0
        if Face_id_counter(Face_id)== 0
            nbrs = neighbors(TR,Face_id);
            nbrs = nbrs(~isnan(nbrs));
            nbrs = (nonzeros(nbrs.*(~Face_id_counter(nbrs))))'; % finds unique faces not yet covered
            Que = [Que nbrs];
            Face_id_counter(Face_id) = 1;
        end
        Face_id = Que(1);
        Que(1) = [];
    end
    
    % Parameters for finding intersection of the circles
    embedded_indices = find(Flags(Face_id,:)==0);
    V1 = Faces(Face_id,embedded_indices(1));
    V2 = Faces(Face_id,embedded_indices(2));
    V3 = nonzeros(Flags(Face_id,:));
    x1 = Flat_V(V1,1);
    y1 = Flat_V(V1,2);
    x2 = Flat_V(V2,1);
    y2 = Flat_V(V2,2);
    r1 = EL(Face_id,embedded_indices(2));
    r2 = EL(Face_id,embedded_indices(1));
    [xout,yout] = circcirc(x1,y1,r1,x2,y2,r2);
   
    %*********
    % DEBUG for NaN's: indicating error in calculating the position of node
    if nnz(isnan([xout;yout]))
        display('error: NaNs in circcirc calculations')
    end
    %*********

    % Choose between points (for when embedding onto a planar mesh)
    % based on the theory that the new point should have the greatest sum
    % of distances from the vertices of previous face (to avoid overlaps or
    % maintain normal vector direction)
    adj_face = cell2mat(edgeAttachments(TR,V1,V2));
    adj_face(adj_face == Face_id) = [];
    pf = Flat_V(Faces(adj_face,:),:)';
    k1 = [xout(1);yout(1)];
    k2 = [xout(2);yout(2)];
    dk1 = sum(sqrt(sum((bsxfun(@minus, pf, k1)).^2)));
    dk2 = sum(sqrt(sum((bsxfun(@minus, pf, k2)).^2)));
    if dk1>dk2
        next_vertex = k1;
    else
        next_vertex = k2;
    end
    
    % Embedding new vertex coordinates into the list
    Flat_V(V3,:) = next_vertex';
    
    % Flagging new vertex as 'embedded'
    Flags(Flags==V3) = 0;
    Face_id_counter(Face_id) = 1;
end