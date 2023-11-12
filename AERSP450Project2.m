clc
clear
close all

%% Analytic Part

% values given in part B analytical part of project promt
a_1 = 13000;
a_2 = 7226.58;              % semi major axes of orbits 1 and 2

e_1 = 0.3;
e_2 = 0.444819;             % eccentricity of orbits 1 and 2

I_1 = 20;
I_2 = 20;                   % inclinations of orbits 1 and 2

RAAN_1 = 30;
RAAN_2 = 30;                % right ascension of the ascending node for orbits 1 and 2

AOP_1 = 50 * pi/180;
AOP_2 = 301.901 * pi/180;   % arguments of periapsis for orbits 1 and 2
                            % given values converted from deg to rad

%calculation of semi latus rectum for both orbits
p_1 = a_1 * (1 - e_1^2);        
p_2 = a_2 * (1 - e_2^2); 

% calculation of gamma, beta, and alpha (given in project guidelines for
% simplified calculations)
gamma = e_1 * e_2 * sin(AOP_1 - AOP_2);
beta = e_1 * p_2 - e_2 * p_1 * cos(AOP_1 - AOP_2);
alpha = e_2 * cos(AOP_1 - AOP_2) - e_1;

% calculation of a, b, and c (determined from equating r^2sin^2f expressions)
a = (e_1^2 - 1) / e_1^2 - alpha^2 / gamma^2;
b = 2*p_1 / e_1^2 - 2*beta*alpha / gamma^2;
c = - (p_1^2 / e_1^2 + beta^2 / gamma^2);


%calculation of radius for orbits 1 and 2 using quadratic formula
r_1 = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
r_2 = (-b - sqrt(b^2 - 4*a*c)) / (2*a)


% calculation of true anaomaly at r - maybe need?
f_1 = acos((((a_1 * (1 - e_1^2)) / r_1) - 1) / e_1);
f_2 = acos((((a_2 * (1 - e_2^2)) / r_2) - 1) / e_2);

%%% from proj 1 idkkk didnt work

% R_P(1) = r_1*cosd(f_1) ;
% R_P(2) = r_1*sind(f_1) ;
% R_P(3) = 0 ;
% 
% cEP = dcm3axis(AOP_1)*dcm1axis(I_1)*dcm3axis(RAAN_1) ;
% cPE = cEP';
% R_ECI = cPE * R_P' 
% % V_ECI = cPE * V_Perifocal' ;
% 
% 
% function r = dcm1axis(ang)
% r = [1 0 0 ; 0 cosd(ang) sind(ang) ; 0 -sind(ang) cosd(ang)];
% end
% 
% % creates a dcm for an angle about axis 3
% function r = dcm3axis(ang)
% r = [cosd(ang) sind(ang) 0 ; -sind(ang) cosd(ang) 0 ; 0 0 1];
% end


%% Numerical Part
load('IODMeasurements.mat')
load('IODMeasurements2.mat')

l = length(AZIMUTH) ;
LOS = zeros(l,3) ;

MU = 3.986*10^5 ;

for i = 1:l
    LOS(i,1) = cosd(ELEVATION(i)) * cosd(AZIMUTH(i)) ;
    LOS(i,2) = cosd(ELEVATION(i)) * sind(AZIMUTH(i)) ;
    LOS(i,3) = sind(ELEVATION(i)) ;
end

RSAT = zeros(3,3,3,l/3) ;
V2SAT = zeros(3,3,l/3) ;
ORBEL = zeros(3,6,l/3) ;
for i = 1:l/3
    [RSAT(:,:,:,i),V2SAT(:,:,i),ORBEL(:,:,i)] = gauss(TIMES((i*3-2:i*3)) ,RSITES((i*3-2:i*3),:) ,LOS((i*3-2:i*3),:),MU) ;
    
end 

init = zeros(1,6) ;
f = figure ;
for i = 1:l/3
    for j = 1:3 
        totalT = 2*pi*sqrt(ORBEL(i,1)^3/MU) ;
        t = linspace(1,totalT,10000) ;
        options = odeset('reltol',1e-12,'abstol',1e-12) ;
        init((1:3)) = RSAT(2,:,:,i) ;
        init((4:6)) = V2SAT(i,:) ;
        [t,RSAT_T] = ode45( @(t,RSAT_T) TwoBP(t,RSAT_T,MU) , t , init, options) ;
        
        
        subplot(1,1,1)
        plot3(RSAT_T(:,1), RSAT_T(:,2), RSAT_T(:,3))
        xlabel('X (KM)')
        ylabel('Y (KM)')
        zlabel('Z (KM)')
        exportgraphics(f,['3D' '.jpg'])
        hold on
    end
end
hold off

function [Rsat,V2sat,OE] = gauss(time,R,L,mu)
    z1 = time(1) - time(2) ;
    z3 = time(3) - time(2) ;
    z13  = time(3) - time(1) ;
    
    lcross(1,:) = cross(L(2,:),L(3,:)) ;
    lcross(2,:) = cross(L(1,:),L(3,:)) ;
    lcross(3,:) = cross(L(1,:),L(2,:)) ;

    D0 = dot(L(1,:),lcross(1,:)) ;
    
    D = zeros(3,3) ;
    for i = 1:3
        for j = 1:3
            D(i,j) = dot(R(i,:),lcross(j,:)) ;
        end
    end

    A = 1/D0 * (-D(1,2)*(z3/z13) + D(2,2) + D(3,2)*(z1/z13)) ;
    B = 1/(6*D0) * ( D(1,2)*(z3^2-z13^2)*(z3/z13) + D(3,2)*(z13^2-z1^2)*(z1/z13)) ;
    E = dot(L(2,:),R(2,:)) ;
    
    R2 = norm(R(2,:)) ;
    a = - (A^2 + 2*A*E + R2) ;
    b = - 2*mu*B*(A+E) ;
    c = - mu^2*B^2 ;

    r2 = roots([1 0 a 0 0 b 0 0 c]) ;
    % rejects values that are not real and positive
    r2 = r2(r2==real(r2)) ;
    r2 = r2(r2>0) ;
    
    Rsat = zeros(3,3,3) ;
    V2sat = zeros(3,3) ;
    OE = zeros(3,6) ;
    for i = 1:length(r2)
        u2 = mu/(r2(i)^3) ;

        c1 = z3/z13 * (1 + (u2/6) * (z13^2 - z3^2)) ;
        c3 = - z1/z13 * (1 + (u2/6) * (z13^2 - z1^2)) ;

        rho(1) = 1/D0 * ( -D(1,1) + (1/c1)*D(2,1) - (c3/c1)*D(3,1) ) ;
        rho(2) = A + u2*B;
        rho(3) = 1/D0 * ( (-c1/c3)*D(1,3) + (1/c3)*D(2,3) - D(3,3) )  ;

        Rsat(1,:,i) = R(1,:) + rho(1) * L(1,:) ;
        Rsat(2,:,i) = R(2,:) + rho(2) * L(2,:) ;
        Rsat(3,:,i) = R(3,:) + rho(3) * L(3,:) 
        
        V2sat(i,:) = gibbs(Rsat(:,:,i),mu) 
       
        [a,e,I,O,W,f] = RV2OE(Rsat(2,:,i),V2sat(i,:),mu) ;
        OE(i,:) = [a,e,I,O,W,f] 
    end
end



function v2 = gibbs(R, mu)
    r(1) = norm(R(1,:)) ;
    r(2) = norm(R(2,:)) ;
    r(3) = norm(R(3,:)) ;
    
    cR2R3 = cross(R(2,:),R(3,:)) ;
    cR3R1 = cross(R(3,:),R(1,:)) ;
    cR1R2 = cross(R(1,:),R(2,:)) ;
       
    D = cR2R3 + cR3R1 + cR1R2 ;
    N = r(1) * cR2R3 +  r(2) * cR3R1 + r(3) * cR1R2 ;
    S = (r(2)-r(3))*R(1,:) + (r(3)-r(1))*R(2,:) + (r(1)-r(2))*R(3,:) ;
    
    d = norm(D) ;
    n = norm(N) ;
    s = norm(S) ;

    v2 = 1/r(2) * sqrt( mu / (n*d)) * cross(D,R(2,:)) + sqrt(mu/(n*d)) * S ;
end



% Two body problem function
function dx = TwoBP(~, r, mu)
    x = r(1) ;
    y = r(2) ;
    z = r(3) ;
    xdot = r(4) ;
    ydot = r(5) ;
    zdot = r(6) ;

    R = sqrt(x^2 + y^2 + z^2) ;

    xdoubledot = -mu/R^3 * x ;
    ydoubledot = -mu/R^3 * y ;
    zdoubledot = -mu/R^3 * z ;

    dx = [ xdot ; ydot ; zdot ; xdoubledot ; ydoubledot ; zdoubledot ] ;
end



function [a,e,I,RAAN,AOP,f] = RV2OE(r,v,mu)
    % declaring gravitational constant of earth and time provided and ECI
    ECI = [[1 0 0];[0 1 0];[0 0 1]];
    
    R1 = norm(r);
    V1 = norm(v);

    energy = (V1^2)/2-mu/R1;
    
    h = cross(r,v);
    H = norm(h);
    
    p = H^2/mu;
    
    % calculating semi major axis
    a = -mu/(2*energy);
    
    % calculating eccentricity 
    eV = cross(v,h)/mu - (r/R1);
    e = norm(eV);
    
    % calculating orbital inclination
    I = acos(dot(ECI(3,:),h)/H);
    
    % calculating longitude of the ascending node
    n = cross(ECI(3,:),h)/norm(cross(ECI(3,:),h));
    
    if(dot(n,ECI(1,:)) >= 0)
        RAAN = atan(dot(n,ECI(2,:))/dot(n,ECI(1,:)));
    elseif (dot(n,ECI(1,:)) < 0)
        RAAN = atan(dot(n,ECI(2,:))/dot(n,ECI(1,:)))+pi;
    end
    
    % calculating argument of periapsis
    if(dot(eV,ECI(:,3)) >= 0)
        AOP = acos(dot(eV,n)/e);
    elseif (dot(eV,ECI(:,3)) < 0)
        AOP = -acos(dot(eV,n)/e);
    end

    if(dot(r,eV) >= 0)
        f = acosd(dot(eV,r)/(e*R1)) ;
    elseif(dot(r,eV) < 0)
        f = 360 - acosd(dot(eV,r)/(e*R1)) ;
    end
end