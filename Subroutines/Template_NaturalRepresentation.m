function [ tri,bary ] = Template_NaturalRepresentation( TRt,DT_FemTemplate )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = DT_FemTemplate.Points;
po = p;
pu = p;
po(:,2) = po(:,2) + 2*pi;
pu(:,2) = pu(:,2) - 2*pi;

%Natural representation
[tri,bary] = pointLocation(TRt,p);

indx = find(isnan(tri));
[t,b] = pointLocation(TRt,po(indx,:));
tri(indx) = t;
bary(indx,:) = b;

indx = find(isnan(tri));
[t,b] = pointLocation(TRt,pu(indx,:));
tri(indx) = t;
bary(indx,:) = b;

% detect NaN nodes that should not overlap
indx = setdiff(indx,find(p(:,1)<min(TRt.Points(:,1))));
%sprintf('There are %d problem nodes: Fixing...',length(indx))
% p1 = plot(p(indx,1),p(indx,2),'.r');

% force NaN nodes that should overlap
er = 1e-6;
noise = [er 0;-er 0;0 er;0 -er];
for i = 1:length(indx)
    pt = p(indx(i),:);
    noise_ou = repmat(noise,[3,1]);
    noise_ou(:,2) = noise_ou(:,2) + reshape((repmat([0,2*pi,-2*pi],[4,1])),12,1);
    temp = bsxfun(@plus,noise_ou,pt);
    [t,b] = pointLocation(TRt,temp);
    if nnz(~isnan(t))
        b(isnan(t),:) = [];
        t(isnan(t)) = [];
    else 
        Problem_node = indx(i)
        %[vi,d] = nearestNeighbor(TRt,[p(indx(i),:);po(3,:);pu(3,:)])
    end
    tri(indx(i)) = t(1);
    bary(indx(i),:) = b(1,:);
end

% Double checks if all NaN's are legitimate (i.e. due to differences in
% shaft lengths and not due to errors in "pointLocation"
final_nan_count = nnz(isnan(tri));
shaft_mismatch_node = nnz(p(:,1)<min(TRt.Points(:,1)));
error = final_nan_count - shaft_mismatch_node;
if error~=0
    message = sprintf('There are still %d problem nodes',length(indx));
    display(message)
else
    %display('No errors to report')
end