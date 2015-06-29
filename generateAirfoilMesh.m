%GENERATE 2D MESH FOR A NACA AIRFOIL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GRID PARAMETERS
L = 0.1; %Length of computational domain (m) 
FineGrid = 512; %Number of points on finest grid
ds = L/(FineGrid*2); %Distance between grid points on finest grid

%AIRFOIL PARAMETERS
first = 2; %NACA NUMBERS (for more info see wikipedia page on naca airfoils)
second = 4;
third = 1;
fourth = 2;
thirdFourth = 12; %last two digits
c = 0.03; %chord length

t = thirdFourth/100.0; %max thickness of the airfoil as a fraction of the chord
m = first/100.0; %maximum camber
p = second/10.0; %location of the maximum camber

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%WRITE .VERTEX FILE
vertex_fid = fopen(['naca2D_' num2str(N) '.vertex'], 'w');

numberNodes = (2*(ceil(c/ds)));

%first line is the number of vertices in the file
fprintf(vertex_fid, '%d\n', numberNodes - 2);
hold on
%remaining lines are the initial coordinates of each vertex

initialX = 0.03;

%WRITE VERTICES FOR UPPER SURFACE
for i = 0:(c/ds)
    %determine x and y coordinates of point along lower surface
    
    in = i*ds; %position along the x axis
    
    
    %HALF THICKNESS OF AIRFOIL AT GIVEN POSITION ON CHORD
    yT = 5*t*c*((0.2969*sqrt(in/c))+(-0.1260*(in/c))+(-0.3516*(in/c)^2)+(0.2843*(in/c)^3)+(-0.1036*(in/c)^4));
    
    %MEAN CAMBER LINE    
    if in <= p*c %calculating for a cambered 4-digit naca airfoil (see wikipedia page)
        
       yC = m*(in/p^2)*(2*p - (in/c));
       
       dyCdx = (2*m/p^2)*(p-(in/c));
       
    else
        
        yC = m*(((c-in)/(1-p)^2)*(1+(in/c)-2*p));
        
        dyCdx = ((2*m/((1-p)^2)))*(p-(in/c));
        
    end
    
    theta = atan(dyCdx);
    
    %UPPER SURFACE OF AIRFOIL
    X(1) = initialX + in - yT*sin(theta); %x coordinates of upper surface
    X(2) = yC + yT*cos(theta); %y coordinates of upper surface
    
    %plot this point
    plot(X(1),X(2),'*r')
    axis([-0.05,0.15,-.05,.05]);
    
    %write the coordinates to the vertex file
    fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
end

%WRITE VERTICES FOR LOWER SURFACE
for i = 1:((c/ds)-1)
    %determine x and y coordinates of point along lower surface
    
    in = i*ds; %position along the x axis
    
    %HALF THICKNESS OF AIRFOIL AT GIVEN POSITION ON CHORD
    yT = 5*t*c*((0.2969*sqrt(in/c))+(-0.1260*(in/c))+(-0.3516*(in/c)^2)+(0.2843*(in/c)^3)+(-0.1036*(in/c)^4));
    
    %MEAN CAMBER LINE    
    if in <= p*c %calculating for a cambered 4-digit naca airfoil (see wikipedia page)
        
       yC = m*(in/p^2)*(2*p - (in/c));
       
       dyCdx = (2*m/p^2)*(p-(in/c));
       
    else
        
        yC = m*(((c-in)/(1-p)^2)*(1+(in/c)-2*p));
        
        dyCdx = ((2*m/((1-p)^2)))*(p-(in/c));
        
    end
    
    theta = atan(dyCdx);
    
    %LOWER SURFACE OF AIRFOIL
    X(1) = initialX + in + yT*sin(theta); %x coordinates of lower surface
    X(2) = yC-yT*cos(theta); %y coordinates of lower surface
    
    %plot this point
    plot(X(1),X(2), '.b')
    
    %write the coordinates to the vertex file
    fprintf(vertex_fid, '%1.16e %1.16e\n', X(1), X(2));
end

hold off
fclose(vertex_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%WRITE .TARGET FILE

targetForce = 1.0e-2;

target_fid = fopen(['naca2D_' num2str(N) '.target'], 'w');

fprintf(target_fid, '%d\n', numberNodes - 2);

for s = 0:numberNodes - 3
   
    fprintf(target_fid, '%d %1.16e\n', s, targetForce);
    
end

fclose(target_fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%