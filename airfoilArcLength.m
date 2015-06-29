%SOME CODE ON FINDING THE ARC LENGTH OF AN AIRFOIL

%NACA NUMBERS (for more info see Wikipedia page on naca airfoils)
first = 2;
second = 4;
third = 1;
fourth = 2;
thirdFourth = 12; %last two digits

c = 0.03; %chord length

ds = c/(N/2); %spatial step for lagrangian points
stepsOverChord = c/ds; %number of steps over chord length

t = thirdFourth/100.0; %max thickness of the airfoil as a fraction of the chord

m = first/100.0; %maximum camber
p = second/10.0; %location of the maximum camber

%INITIALIZE SOME VARIABLES

%COORDINATES OF LEADING TIP
x0 = 0.0;
y0 = 0.0;

%ARC LENGTH
arcLengthUpper = 0; %arc length of upper surface
arcLengthLower = 0; %arc length of lower surface
arcLengthTotal = 0; %arc length of airfoil surface

%SET FIRST SET OF "PREVIOUS" COORDINATES TO ORIGIN OF AIRFOIL
xUPrevious = 0.0;
yUPrevious = 0.0;

xLPrevious = 0.0;
yLPrevious = 0.0;

%%CREATE UPPER AND LOWER SURFACES OF AIRFOIL...
for x = 1:N
    
    xn = x*dx; %position along the x axis
    
    %HALF THICKNESS OF AIRFOIL AT GIVEN POSITION ON CHORD
    yT = 5*t*c*((0.2969*sqrt(xn/c))+(-0.1260*(xn/c))+(-0.3516*(xn/c)^2)+(0.2843*(xn/c)^3)+(-0.1015*(xn/c)^4));
    
    %MEAN CAMBER LINE    
    if xn <= p*c %calculating for a cambered 4-digit naca airfoil (see wikipedia page)
        
       yC = m*(xn/p^2)*(2*p - (xn/c));
       
       dyCdx = (2*m/p^2)*(p-(xn/c));
       
    else
        
        yC = m*(((c-xn)/(1-p)^2)*(1+(xn/c)-2*p));
        
        dyCdx = ((2*m/((1-p)^2)))*(p-(xn/c));
        
    end
    
    theta = atan(dyCdx);
    
    %UPPER SURFACE OF AIRFOIL
    yU(x) = yC + yT*cos(theta); %y coordinates of upper surface
    xU(x) = xn - yT*sin(theta); %x coordinates of upper surface
    
    %find distance from this point to previous
    dsUpper = sqrt((xU(x) - xUPrevious)^2 + (yU(x) - yUPrevious)^2);
    
    %update future previous values (x(i-1))
    xUPrevious = xU(x);
    yUPrevious = yU(x);

    %add the step in length to the total arc length (upper)
    arcLengthUpper = arcLengthUpper + dsUpper;
    
    %LOWER SURFACE OF AIRFOIL
    xL(x) = xn + yT*sin(theta); %x coordinates of lower surface
    yL(x) = yC-yT*cos(theta); %y coordinates of lower surface
    
    %find distance from this point to previous
    dsLower = sqrt((xL(x) - xLPrevious)^2 + (yL(x) - yLPrevious)^2);
    
    %update future previous values (x(i-1))
    xLPrevious = xL(x);
    yLPrevious = yL(x);
    
    %add the step in length to the total arc length (lower)
    arcLengthLower = arcLengthLower + dsLower;
    
end

%find total arc length of the airfoil
arcLengthTotal = arcLengthLower + arcLengthUpper;

%PLOT THE GEOMETRY
%plot(xU, yU, 'r-'); hold on; %upper surface of airfoil
%plot(xU, yU, '*'); hold on;

%plot(xL, yL, 'r-'); hold on; %lower surface of airfoil
%plot(xL, yL, '*'); hold on;

%xlabel('x'); ylabel('y');
%axis([-0.05,0.15,-.05,.05]);

%FIND THE NUMBER OF NODES TO LAY ALONG AIRFOIL
numNodesUpper = ceil(arcLengthUpper/(0.5*dx)/4)*4+1;
numNodesLower = ceil(arcLengthLower/(0.5*dx)/4)*4+1;
numNodesTotal = numNodesUpper + numNodesLower;

%FIND LENGTH STEP ALONG ARC
dsUpper = arcLengthUpper/(numNodesUpper-1);
dsLower = arcLengthLower/(numNodesLower-1);

kappa_target = 1.0e-2; %target point penalty spring constant (Newton)
