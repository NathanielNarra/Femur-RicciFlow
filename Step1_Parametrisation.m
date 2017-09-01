function [] = Step1_Parametrisation( file_paths )
% Read-in surface meshes one by one and parametrise as an annulus.
%   Detailed explanation goes here
% INPUT:'file_paths' containing the full path to the list of surface meshes
%   examples-
%   [~,file_paths]=fileattrib('*.obj');
%   [~,file_paths]=fileattrib('*.ply');

for f = 1:length(file_paths)
    [V,F] = read_ply(file_paths(f).Name);
    TR = triangulation(F,V);
    NofVertices = size(TR.Points,1);
    
    % Free boundary conformal parametrization
    [EL,ss] = Discrete_Ricci_Flow(TR,1);
    flat = Seed_face(EL,TR);
    TR_flat = triangulation(TR.ConnectivityList,flat);
    cfactor = ss.U0- ss.U(:,end);
    cfactor_sm = cfactor;
    
    % smooth conformal factor 'cfactor' over a 2 ring neighbourhood
    R = Ring_Neighbours(TR,2);
    for i = 1:NofVertices
        cfactor_sm(i) = mean(cfactor(cell2mat(R(:,i))));
    end
    
    % Feature points: Femoral Head (FH) and Greater Trochanter (GT)
    [FHnode,GTnode,FHb,FHsupport] = Detect_Landmarks(TR_flat,cfactor_sm);
    
    % Manipulates the mesh so nodes lie on shortest path FH-GT 
    [TR_flat,TR,v_ix,f_ix] = Path_straighten(TR_flat,TR,FHnode,GTnode);
    cfactor_sm = cfactor_sm(v_ix);
    NofVertices = size(TR.Points,1);
    orig_bound = unique(freeBoundary(TR));
    FHnode = f_ix(FHnode);
    GTnode = f_ix(GTnode);
    FHb = f_ix(FHb);
    FHsupport = f_ix(FHsupport);
    
    % Shortest path between the 2 feature points
    E = edges(TR);
    NofEdges = size(E,1);
    VconnE = sparse([E(:,2);E(:,1)],[E(:,1);E(:,2)],[(1:NofEdges)';(1:NofEdges)'],NofVertices,NofVertices);
    
    % Shortest FH-GT path in disk parametrisation
    [~, Edge_Lengthsf] = Discrete_Riemannian_Metric(TR_flat);
    G = spfun(@(x) Edge_Lengthsf(x),VconnE);
    [~, slitpath, ~] = graphshortestpath(G,FHnode,GTnode);
    
    % Check that the detected feature nodes are not of valence:3. This
    % creates problems in the slit creation if the terminal nodes share the
    % same number of neighbours as the next node in path
    Valence = sum(spones(VconnE));
    if Valence(FHnode)==3
        FHnode = slitpath(2); 
        slitpath(1) = [];
        disp('FHnode detected to be Valence3. Path corrected accordingly')
    end
    if Valence(GTnode)==3
        GTnode = slitpath(end-1);
        slitpath(end) = [];
        disp('GTnode detected to be Valence3. Path corrected accordingly')
    end
    
    % Slit mesh topology along the shortest path FH-GT
    [Fn,Vn] = MeshSlit(TR,slitpath',1);
    TR_double = triangulation(Fn,Vn);
    
    % re-parametrize the new doubly connected Genus 0 mesh. Target
    % curvatures set to 0 every where (including the 2 boundaries)
    [ELn,ssn] = Discrete_Ricci_Flow(TR_double,0);
    cf = ssn.U0- ssn.U(:,end);
    
    % find shortest path from GTnode to original boundary wrt to original
    % mesh embedding & slit the mesh along this path
    [dist, cutgraph_path, ~] = graphshortestpath(G,GTnode,orig_bound);
    
    % check that the path between the GT feature point and the shaft
    % boundary does not share a path with the inter-feature slitpath
    j = 1;
    [~,sorted_dist_indx] = sort(dist);
    sp = cutgraph_path{sorted_dist_indx(1)};
    while ~isempty(intersect(sp(2:end-1),slitpath(2:end-1)))
        j = j+1;
        sp = cutgraph_path{sorted_dist_indx(j)};
    end
    [F_slit,V_slit] = MeshSlit(TR_double,sp');
    TR_slit = triangulation(F_slit,V_slit);
        
    % embed the mesh in the parametric domain and update conformal factor
    % list to account for additional nodes (due to cutgraph)
    flat_slit = Seed_face(ELn,TR_slit);
    cdata_slit = [cf; cf(sp)];
    % Detect the inter-feature boundary (FH-GT)
    [boundary1,boundary2] = Detect_Boundaries(TR_double);
    if isequal(sort(boundary1),sort(orig_bound))
        B = boundary2;
    else
        B = boundary1;
    end
    an = setdiff(find(ismember(TR_slit.Points,TR_slit.Points(B,:),'rows')),B);
    B = [B;an];
    [Annulus,Z4] = Annulus_ComplexPlane(flat_slit,B);
    
    % Update FHboundary list to include duplicated nodes
    FHboundary = find(ismember(round(V_slit*10000)/10000,round(TR.Points(FHb,:)*10000)/10000,'rows'));
    FHsupport = find(ismember(round(V_slit*10000)/10000,round(TR.Points(FHsupport,:)*10000)/10000,'rows'));
    
    % self referencing file naming system (change, if not ply file)
    basename = strrep(file_paths(f).Name,'.ply','');
    filename = sprintf('%s_Mapped.mat',basename);
    
    % Save data structure for each sample onto disk
    Parametrization.Simple = TR;
    Parametrization.Disk = TR_flat;
    Parametrization.CFdisk = cfactor_sm;
    Parametrization.Double = TR_double;
    Parametrization.DoubleSlit = TR_slit;
    Parametrization.FHnode = FHnode;
    Parametrization.FHboundary = FHboundary;
    Parametrization.FHsupport = FHsupport;
    Parametrization.GTnode = GTnode;
    Parametrization.ParamPlane = Z4;
    Parametrization.Annulus = Annulus;
    Parametrization.ConFac = cdata_slit;
    save(filename,'Parametrization');
    
    % OPTIONAL: Display the FH support regions and the inter-feature
    % boundary to monitor feature detection related processes
    figure, trimesh(TR_slit); axis equal;set(gca,'Zdir','reverse'), hold on
    plot3(V_slit(FHnode,1),V_slit(FHnode,2),V_slit(FHnode,3),'r.','MarkerSize',20);
    plot3(V_slit(FHboundary,1),V_slit(FHboundary,2),V_slit(FHboundary,3),'g.','MarkerSize',20);
    plot3(V_slit(FHsupport,1),V_slit(FHsupport,2),V_slit(FHsupport,3),'c.');
    plot3(V_slit(B,1),V_slit(B,2),V_slit(B,3),'k.');
    hold off
    savefig(basename)
    close
    
    disp(sprintf('Subject %d mesh mapped and data saved',f))
end
end