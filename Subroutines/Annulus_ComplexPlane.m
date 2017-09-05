function [ Annulus,Z3 ] = Annulus_ComplexPlane(flat,B)
%Normalise the parametrised mesh in the coordinate plane
%   B: outer boundary nodes 

Z = complex(flat(:,1),flat(:,2));
Zb = Z(B);
[~,ix] = min(real(Zb));
Z1 = Z - Zb(ix);
%figure, f = trimesh(TR.ConnectivityList,real(Z1),imag(Z1),cdata_slit), axis equal
Z1b = Z1(B);
theta = median(angle(Z1b));
Z2 = Z1*exp(1i*(pi/2-theta));

%scaling:2pi length
sf = max(imag(Z2(B)));
Z3 = ((2*pi)/sf)*Z2;

% Place consistent boundary (B) on the imaginary axis
if mean(real(Z3)) > 0 
    Z3 = -1*conj(Z3); % flip about the imaginary axis
end

if max(abs(real(Z3(B)))) > 1e-6
    warning('Boundary nodes not sufficiently aligned with imaginary axis, ERROR: %d',max(abs(real(Z3(B)))));
end
% enforcing real component of boundary nodes B to be absolute zero
Z3(B) = (Z3(B) - conj(Z3(B)))/2;
%Planar annulus
Annulus = exp(Z3);
end

