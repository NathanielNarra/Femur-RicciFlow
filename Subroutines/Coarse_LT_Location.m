function [ lt_target, DT_target ] = Coarse_LT_Location( roi,lt_template, DT_target,cf_target )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

X = DT_target.Points(:,1);
Y = DT_target.Points(:,2);
Z = cf_target;

% wrap nodes outside [0,2*pi]
Y(Y>2*pi) = Y(Y>2*pi) - 2*pi;
Y(Y<0) = Y(Y<0) + 2*pi;

Fn = scatteredInterpolant(X,Y,Z,'linear','nearest');
[gridx, gridy] = meshgrid(-(2:pi/100:max(abs(X))),min(Y):pi/100:max(Y));
gridx = flip(gridx,2);
gridy = flip(gridy,1);
grid = Fn(gridx,gridy);
grid = grid - mean(grid(:));

corr_img = xcorr2(grid,roi);
[r,c] = size(roi);
[rr,cc] = size(grid);
ro = floor(r/2);
co = floor(c/2);
I = corr_img(ro+1:ro+rr,co+1:co+cc);
I_flipped = flip(I,1);

[row,col] = find(I==max(I(:)));
[rowf,colf] = find(I_flipped==max(I_flipped(:)));

lt_target = [gridx(row,col), gridy(row,col)];
lt_target_f = [gridx(rowf,colf), gridy(rowf,colf)];

if pdist2(lt_target,lt_template) > pdist2(lt_target_f,lt_template)
    lt_target = lt_target_f;
    DT_target = triangulation(DT_target.ConnectivityList,DT_target.Points(:,1),2*pi-DT_target.Points(:,2));
end
end

