 function [EL_final, SS] = Discrete_Ricci_Flow(TR,FB)

% Calculate the Riemannian metric
[RM, Edge_Lengths] = Discrete_Riemannian_Metric(TR);

% Calculate the initial circle metric
[I, Gamma_f, Gamma_v,I_list]= Initial_Circle_Packing_Metric(TR,RM,Edge_Lengths);

% Evolve the metric
[EL_final, SS] = Newton_Method(TR,I,Gamma_f,Gamma_v,I_list,FB);

end


function [Edge_Weights, Gamma_f, Gamma_v, InvD] = Initial_Circle_Packing_Metric(TR,RM,EL)
TRI = TR.ConnectivityList;

% Compute Gamma values for each vertex for each Triangle (Face)
Gamma_f = (repmat(sum(RM,2),[1,3])-2*RM)/2;

% Compute Gamma values for each vertex as a function of its multiple values
% from adjacent triangles
Gamma_v = struct2array(regionprops(TRI,Gamma_f,'MinIntensity'))';

% Compute inversive distance
temp = Gamma_v(edges(TR));
InvD = (EL.^2 - sum(temp.^2,2))./(2*prod(temp,2));
Gamma_f = Gamma_v(TRI);
Gi = circshift(Gamma_f,[0 -1]);
Gj = circshift(Gamma_f,[0 -2]); 
Edge_Weights = (RM.^2 - Gi.^2 - Gj.^2)./(2*Gi.*Gj);
end

function [EL, SS] = Newton_Method(TR,I,~,Gamma_v,I_list,FreeBoundary)
noV = size(TR.Points,1);
U = log(Gamma_v);
E = edges(TR);
SS = {};
SS.U0 = U;

Gamma_v = exp(U);
Gamma_f = Gamma_v(TR.ConnectivityList);


%**** Potential to optimse 
FE = Faces_constituting_Edges(TR);
Edge_indices_in_F = Edge_Instances(TR,FE);
[V_matrix_edgeid,~,V_Occur,~] = Vertex_Connectivity(TR);
%****

b_vertices = unique(freeBoundary(TR));
interior_vertices = setdiff(unique(TR.ConnectivityList),b_vertices);

% Current curvature
Gi = circshift(Gamma_f,[0 -1]);
Gj = circshift(Gamma_f,[0 -2]);
EL = sqrt(Gi.^2 + Gj.^2 + 2*(Gi.*Gj.*I));
numerator = repmat(sum(EL.^2,2),[1,3]) - 2*(EL.^2);
denominator = repmat(2*prod(EL,2),[1,3])./EL;
corner_angles =  acos(numerator./denominator);
Theta_sums = full(sum(spfun(@(x) corner_angles(x),V_Occur))');
K = 2*pi - Theta_sums;
K(b_vertices) = pi - Theta_sums(b_vertices);

% Target Curvature
K_target = zeros(size(Gamma_v)); 
Error = max(abs(-K(interior_vertices))); % if only internal vertex metrics are changed
iteration_count = 0;
while Error > 0.000001
    
    % Edge lengths
    Gi = circshift(Gamma_f,[0 -1]);
    Gj = circshift(Gamma_f,[0 -2]);
    EL = sqrt(Gi.^2 + Gj.^2 + 2*(Gi.*Gj.*I));
    
    % Computing corner angles
    numerator = repmat(sum(EL.^2,2),[1,3]) - 2*(EL.^2); 
    denominator = repmat(2*prod(EL,2),[1,3])./EL;
    Q = numerator./denominator;
    % diagnose where the triangles are degenerating
    if max(Q) >1 | min(Q)<-1
        trimesh(TR), axis equal
        hold on
        ix = find(any(Q>1,2));
        p1 = plot3(TR.Points(TR.ConnectivityList(ix,:),1),TR.Points(TR.ConnectivityList(ix,:),2),TR.Points(TR.ConnectivityList(ix,:),3),'*r');
        error('acos out of bounds: probably bad triangle quality')
    end
    
    corner_angles =  acos(numerator./denominator);
    Theta_sums = full(sum(spfun(@(x) corner_angles(x),V_Occur))');
    K = 2*pi - Theta_sums;
    K(b_vertices) = pi - Theta_sums(b_vertices);
    
    % Orthocentre and distances therefrom 
    Orth_dist = Orthocentre_Edge_Distance(TR,EL,corner_angles,Gamma_f);
    
    %Edge weights and Hessian Matrix construction.
    W_matrix = Orth_dist./EL;
    W_edge_weights = full(sum(spfun(@(x) W_matrix(x),Edge_indices_in_F))');
    H = spfun(@(x) W_edge_weights(x),V_matrix_edgeid);
    H = diag(sum(H,2)) - H;

    % For free boundary conditions: truncating the Hessian by ignoring
    % the boundary node entries
    if FreeBoundary
        H(b_vertices,:) = [];
        H(:,b_vertices) = [];
        L = ichol(H); %'H' must be sparse
        dU_int = pcg(H,(K_target(interior_vertices)-K(interior_vertices)),[],noV,L,L');
        U(interior_vertices) = U(interior_vertices) + dU_int;
    else
        % Preconditioning and Conjugate gradient method
        L = ichol(H);
        dU = pcg(H,(K_target - K),[],noV,L,L');
        U = U + dU; % all node metrics are updated
        U = U - mean(U);
    end
    
    Gamma_v = exp(U);
    Gamma_f = Gamma_v(TR.ConnectivityList);
    
    Error = max(abs(-K(interior_vertices))); 
    iteration_count = iteration_count + 1;
    
    % Collecting all data structures
    SS.U(:,iteration_count) = U;
    Gij = Gamma_v(E);
    SS.L(:,iteration_count) = sqrt(sum(Gij.^2,2) + (prod(Gij,2).*I_list));
    SS.K(:,iteration_count) = K;
    SS.CornerAngles(:,:,iteration_count) = corner_angles;    
end
% Final edge lengths
Gi = circshift(Gamma_f,[0 -1]);
Gj = circshift(Gamma_f,[0 -2]);
EL = sqrt(Gi.^2 + Gj.^2 + 2*(Gi.*Gj.*I));
end