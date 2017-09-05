function ortho_distance = Orthocentre_Edge_Distance(TR,EL,angles,G_f)

ortho_distance = zeros(size(EL));
for i = 1:size(TR.ConnectivityList,1)
    A = [EL(i,3) 0;EL(i,2)*cos(angles(i,1)) EL(i,2)*sin(angles(i,1))];
    B = [EL(i,3)^2 + G_f(i,1)^2 - G_f(i,2)^2; EL(i,2)^2 + G_f(i,1)^2 - G_f(i,3)^2];
    oc = 0.5*(A\B);
    
    % create local coordinates (planar,2d) for the 3 vertices of the face
    % size: 2x3. rows -> x and y. columns -> vertices 1,2 & 3
    vertex_coordinates = [0 EL(i,3) EL(i,2)*cos(angles(i,1));0 0 EL(i,2)*sin(angles(i,1))];
    
    t = [oc oc oc] - vertex_coordinates;
    a = sqrt(sum(t.^2,1));
    b = [a(2) a(3) a(1)];
    c = [EL(i,3) EL(i,1) EL(i,2)];
    s = (a+b+c)/2; % 1x3 vector
    h = 2*sqrt(s.*(s-a).*(s-b).*(s-c))./c;
    % in case the ortho centre lies on one of the sides of the triangle,
    % the loop truncates complex numbers to their real parts
    if(~isreal(h)) 
        h = real(h);
        warning('Orthocentre of face lies on one of its sides')
        face = i
    end
    ortho_distance(i,:) = [h(2) h(3) h(1)];
    
end