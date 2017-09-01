function [NatRep,node_coordinates3D, thk] = STEP2_Template_Matching( folder_id )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% % For thickness data
% load('C:\Users\narragir\Documents\MATLAB\Ricci Flow\FemurTemplate.mat')
% cf = Parametrization.ConFac;
% clear Parametrization
% dir_path = fullfile('C:\Users\narragir\Documents\MATLAB\Ricci Flow\Paper4 Data Collection');
% folder = sprintf('Group%d',folder_id);

% For Von Mises data
load('C:\Users\narragir\Documents\MATLAB\Ricci Flow\Paper4 Data Collection\ForNateSurfaceVM\VM_Template.mat')
dir_path = fullfile('C:\Users\narragir\Documents\MATLAB\Ricci Flow\Paper4 Data Collection\ForNateSurfaceVM');
folder = sprintf('%d00',folder_id);
dir_path = fullfile(dir_path,folder);
files = dir(strcat(dir_path,'\vm*Mapped*.mat'));
clear Parametrization

DT_template = triangulation(F,V);
lmloc_template = V(LTnode,:);
%***** Template LT process
d = pdist2(V,V(LTnode,:));
roi = d<0.7;
ltnx = V(roi,1);
ltny = V(roi,2);
ltnz = cf(roi);
[tgridx, tgridy] = meshgrid(min(ltnx):pi/100:max(ltnx),min(ltny):pi/100:max(ltny));
tFn = scatteredInterpolant(ltnx,ltny,ltnz);
LTroi_template = tFn(tgridx,tgridy);
LTroi_template = LTroi_template - mean(LTroi_template(:));
%*****

Subjectcount = [19,19,17,18,18,20];
NofSub = Subjectcount(folder_id);
GTnodes_template = find(ismember(V,[0,0;0,2*pi],'rows'));
bnodes_template = unique(freeBoundary(DT_template));


for subject_id = 1:NofSub
    %file_id
    
    filename_mapped = fullfile(dir_path,files(subject_id).name);
    display(strcat('loading____',files(subject_id).name))
    load(filename_mapped)
    
    cf_target = Parametrization.ConFac;
    DT_target = triangulation(Parametrization.DoubleSlit.ConnectivityList,real(Parametrization.ParamPlane),imag(Parametrization.ParamPlane));
    FHnode_target = Parametrization.FHnode;
    % smooth cf_target: 2 ring neighbours
    ring = Ring_Neighbours(DT_target,2);
    cf_target = cellfun(@(x) mean(cf_target(x)),ring(3,:))';
    
    % ascertain proper orientation of the maps wrt template
    [lmloc_target_rough,DT_target] = Coarse_LT_Location(LTroi_template,lmloc_template,DT_target,cf_target);
    
    % Fitting gaussian surface to LT process
    % In case the LTlocation is split across
    P = DT_target.Points;
    PP = [P;bsxfun(@plus,P,[0,2*pi])];
    
    P(P(:,2)>2*pi,2) = P(P(:,2)>2*pi,2) - 2*pi;
    P(P(:,2)<0,2) = P(P(:,2)<0,2) + 2*pi;
    Vcentered = bsxfun(@minus,P,lmloc_target_rough);
    
    dist = sqrt(sum(Vcentered.^2,2));
    indx = find(dist<1);
    
    weights = (1-dist(indx)).^2;
    x = Vcentered(indx,1);
    y = Vcentered(indx,2);
    cfmod = [cf_target;cf_target];
    z = cfmod(indx);
    z = z - min(z);
    g = fittype('a*exp(-((x-xm)^2/(2*sx^2) + (y-ym)^2/(2*sy^2)))','dependent','z','independent',{'x','y'},'coefficients',{'a','sx','sy','xm','ym'});
    fopt = fitoptions(g);
    fopt.Weights = weights;
    fopt.StartPoint = [0.9*max(z),0.25,0.25,0,0];
    [gaussfit, ~]= fit([x(:),y(:)],z(:),g,fopt);
    
    %     figure, plot(gaussfit,[x,y],z)
    %     figure('Name',['Subject ',num2str(subject_id)]),
    %     %f1 = subplot(1,2,1); plot(fit1,[x,y],z)
    %     f2 = subplot(1,2,2); plot(fit2,[x,y],z)
    %     %title(f1,['Unweighted - rsquare: ',num2str(gof1.rsquare)])
    %     title(f2,['Weighted - rsquare: ',num2str(gof2.rsquare)])
    %     pause
    %     close
    
    lmloc_target = lmloc_target_rough + [gaussfit.xm, gaussfit.ym];
    
    % Match LT region
    % (includes boundary nodes as fixed nodes such that
    % movement of internal nodes does not alter the node positions on the
    % perimeter of the template mesh)
    lm_displacement = lmloc_target - lmloc_template;
    if sqrt(sum(lm_displacement.^2)) > 1.2
        figure, plot(gaussfit,[x,y],z)
        display('surface fitting failure...')
    end
    
    
    % Based on finding new positions for the boundary nodes first
    x_shaft = min(V(:,1));
    scale = sqrt(sum(lm_displacement.^2,2)) * 4;
    shaft_disp = 0.99*lm_displacement(1);
    d = abs(V(bnodes_cut,1) - x_shaft);
    r = d/scale;
    k = ((1-r).^4).*(4*r+1);
    k(r>1) = 0;
    cut_disp = k*shaft_disp;
    
    np_lt = bsxfun(@plus,V(LTsupport,:),lm_displacement);
    np_shaft = bsxfun(@plus,V(bnodes_shaft,:),[shaft_disp,0]);
    np_cut = bsxfun(@plus,V(bnodes_cut,:),[cut_disp,zeros(size(cut_disp))]);
    
    source = [V(LTsupport,:); V(bnodes_shaft,:); V(bnodes_cut,:); V(bnodes_slit,:)];
    target = [np_lt; np_shaft; np_cut; V(bnodes_slit,:)];
    [modified_pos,~] = RBF(V,source,target,scale);
    
    %****************************
    % Match FH region
    fhb_target = Parametrization.FHboundary;
    V_target = DT_target.Points;
    fh_displacement = DT_target.Points(FHnode_target,:) - modified_pos(FHnode,:);
    scale = 4*sqrt(sum(fh_displacement.^2,2));
    source = modified_pos(FHsupport,:);
    target = bsxfun(@plus,source,fh_displacement);
    if scale>0
        [modified_pos,~] = RBF(modified_pos,source,target,scale);
    end
    [ix,~] = knnsearch(V_target(fhb_target,:),modified_pos(FHboundary,:));
    % simple scaling of FHsupport
    p = complex(modified_pos(FHboundary,1),modified_pos(FHboundary,2));
    o = complex(modified_pos(FHnode,1),modified_pos(FHnode,2));
    vectors = p-o;
    direction  = abs(complex(V_target(fhb_target(ix),1),V_target(fhb_target(ix),2)) - o) - abs(vectors);
    mag_factor = mean((abs(vectors)+direction)./abs(vectors));
    % RBF with FHsupport new positions as target
    ps = complex(modified_pos(FHsupport,1),modified_pos(FHsupport,2));
    vs = ps - o;
    result_s = mag_factor*vs + o;
    source = modified_pos(FHsupport,:);
    target = [real(result_s),imag(result_s)];
    scale = 4 * max(abs(ps) - abs(result_s));
    if scale>0
        [modified_pos,~] = RBF(modified_pos,source,target,scale);
    end
    %***************************
    
    modified_pos(GTnodes_template,:) = [0,0;0,2*pi];
    DT_template_mod = triangulation(DT_template.ConnectivityList,modified_pos(:,1),modified_pos(:,2));
    [Y,path_all] = Matching_Slitnodes(DT_template_mod,DT_target,FHnode,FHnode_target);
    
    %RBF
    source = modified_pos(path_all,:);
    target = [zeros(size(Y)),Y];
    max_disp = max(abs(Y-source(:,2)));
    scale = 4*max_disp;
    [modified_pos,~] = RBF(modified_pos,source,target,scale);
    
    % Ensure slit nodes remain on y axis (x=0)
    modified_pos(bnodes_slit,1) = 0;
    modified_pos(GTnodes_template,:) = [0,0;0,2*pi];
    % Ensure slit nodes remain on y axis (x=0)
    modified_pos(bnodes_cut_bottom,2) = 0;
    modified_pos(bnodes_cut_top,2) = 2*pi;
        
    [tri_reg,bary_reg] = Template_NaturalRepresentation(DT_target,triangulation(DT_template.ConnectivityList,modified_pos(:,1),modified_pos(:,2)));
    nan_indx = isnan(tri_reg);
    tri_reg(nan_indx) = 1; %dummy value
    
    %***********************************
    %     % Calculating thickness features
    %     filename_thickness = fullfile(dir_path,folder,sprintf('G%d_Thickness_%d.mat',folder_id,subject_id));
    %     [thickness, thickness_minimum] = Thickness_finetuning( filename_thickness );
    %     temp = thickness(Parametrization.Simple.ConnectivityList(tri_reg,:));
    %     temp1 = thickness_minimum(Parametrization.Simple.ConnectivityList(tri_reg,:));
    %     temp(nan_indx,:) = NaN;
    %     temp1(nan_indx,:) = NaN;
    %     thk.surfnormal(:,subject_id) = sum(temp.*bary_reg,2);
    %     thk.minimum(:,subject_id) = sum(temp1.*bary_reg,2);
    %***********************************
    thk = 'null'; % temporary filler
    
    %temp_position = zeros(size(V3d)); %for thickness analysis
    temp_position = zeros(size(V,1),3);
    temp_position(~nan_indx,:) = barycentricToCartesian(Parametrization.Simple,tri_reg(~nan_indx),bary_reg(~nan_indx,:));
    node_coordinates3D(:,:,subject_id) = temp_position;
    node_coordinates3D(nan_indx,:,subject_id) = NaN;
    NatRep(subject_id).faceID = tri_reg;
    NatRep(subject_id).vertexID = Parametrization.Simple.ConnectivityList(tri_reg,:);
    NatRep(subject_id).barycent = bary_reg;
    NatRep(subject_id).nan_indx = nan_indx;
    NatRep(subject_id).Mod_Template = triangulation(F,modified_pos);
    %NatRep.V3D = node_coordinates3D;
    
    % Verification: Graphs and Plots
    %     lm1 = plot3(f1,lmloc_template(1),lmloc_template(2),1,'.r','MarkerSize',16);
    %     lm2 = plot3(f2,lmloc_template(1),lmloc_template(2),1,'.r','MarkerSize',16);
    %     lm3 = plot3(f2,lmloc_target(1),lmloc_target(2),1,'.b','MarkerSize',16);
    %     PC = barycentricToCartesian(Parametrization.Simple,tri(~isnan(tri)),bary(~isnan(bary)));
    %     figure,
    %     f1 = subplot(1,2,1); trimesh(DT_template.ConnectivityList,DT_template.Points(:,1),DT_template.Points(:,2),cf), axis equal
    %     set(f1,'View',[0 90]);
    %     hold(f1)
    %     f2 = subplot(1,2,2); trimesh(DT_template.ConnectivityList,modified_pos(:,1),modified_pos(:,2),cf), axis equal
    %     hold(f2)
    %     trimesh(DT_target.ConnectivityList,-DT_target.Points(:,1),DT_target.Points(:,2),cf_target), axis equal
    %     set(f2,'View',[0 90]);
    
    figure
    f1 = subplot(1,2,1); trimesh(DT_template.ConnectivityList,node_coordinates3D(:,1,subject_id),node_coordinates3D(:,2,subject_id),node_coordinates3D(:,3,subject_id),cf); axis equal, set(gca,'Zdir','reverse')
    f2 = subplot(1,2,2); P = trimesh(Parametrization.DoubleSlit); axis equal, set(gca,'Zdir','reverse')
    set(f1,'CameraPosition',[702.4480  341.1933 -121.9889]);
    set(f2,'CameraPosition',[702.4480  341.1933 -121.9889]);
    set(P,'FaceColor','interp','FaceVertexCData',cf_target,'EdgeColor',[0 0 0],'EdgeAlpha',0.2,'LineWidth',1,'CDataMapping','scaled')
    hold(f1)
    p1 = plot3(f1,node_coordinates3D(LTsupport,1,subject_id),node_coordinates3D(LTsupport,2,subject_id),node_coordinates3D(LTsupport,3,subject_id),'.r','MarkerSize',10);
    p2 = plot3(f1,node_coordinates3D(FHsupport,1,subject_id),node_coordinates3D(FHsupport,2,subject_id),node_coordinates3D(FHsupport,3,subject_id),'.k','MarkerSize',10);
    fn = sprintf('Group %d Subject %d.fig',folder_id,subject_id);
    title(fn)
    savefig(fn)
    close
    
    %         figure,
    %     f1 = subplot(1,2,1); trimesh(DT_template.ConnectivityList,DT_template.Points(:,1),DT_template.Points(:,2),cf), axis equal
    %     set(f1,'View',[0 90]);
    %     f2 = subplot(1,2,2); trimesh(DT_target.ConnectivityList,DT_target.Points(:,1),DT_target.Points(:,2),cf_target), axis equal
    %     set(f2,'View',[0 90]);
    %     pause
    %     hold(f1);hold(f2);
end

