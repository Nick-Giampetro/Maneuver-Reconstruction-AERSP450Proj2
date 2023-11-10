clc
clear
close all

load('IODMeasurements.mat')
load('IODMeasurements2.mat')

l = length(AZIMUTH) ;
LOS = zeros(l,3) ;

MU = 3.986*10^5

for i = 1:l
    LOS(i,1) = cosd(ELEVATION(i)) * cosd(AZIMUTH(i)) ;
    LOS(i,2) = cosd(ELEVATION(i)) * sind(AZIMUTH(i)) ;
    LOS(i,3) = sind(ELEVATION(i)) ;
end

RSAT = zeros(3,3,l/3) ;
V2SAT = zeros(l/3,3) ;
ORBEL = zeros(l/3,6) ;
for i = 1:l/3
    [RSAT(:,:,i),V2SAT(i,:),ORBEL(i,:)] = gauss(TIMES((i*3-2:i*3)) ,RSITES((i*3-2:i*3),:) ,LOS((i*3-2:i*3),:),MU,Re)
    
end 



function [Rsat,V2sat,OE] = gauss(time,R,Phat,mu,Re)
    t1 = time(1) - time(2) ;
    t3 = time(3) - time(2) ;
    t  = time(3) - time(1) ;
    
    p(1,:) = cross(Phat(2,:),Phat(3,:)) ;
    p(2,:) = cross(Phat(1,:),Phat(3,:)) ;
    p(3,:) = cross(Phat(1,:),Phat(2,:)) ;

    D0 = dot(Phat(1,:),p(1,:)) ;
    
    D = zeros(3,3) ;
    for i = 1:3
        for j = 1:3
            D(i,j) = dot(R(i,:),p(j,:)) ;
        end
    end

    A = 1/D0 * (-D(1,2)*t3/t + D(2,2) + D(3,2)*t1/t) ;
    B = 1/(6*D0) * ( D(1,2)*(t3^2-t^2)*t3/t + D(3,2)*(t^2-t1^2)*t1/t) ;
    E = dot(R(2,:),Phat(2,:)) ;
    
    R2 = dot(R(2,:),R(2,:)) ;
    a = - (A^2 + 2*A*E + R2) ;
    b = - 2*mu*B*(A+E) ;
    c = - mu^2*B^2 ;

    r2 = roots([1 0 a 0 0 b 0 0 c]) ;
    % rejects values that are not positive 
    r2 = r2(r2==real(r2)) ;
    r2 = r2(r2>0) ;
    
    for i = 1:length(r2)
        rho(1) = 1/D0 * ( (6*(D(3,1)*t1/t3 + D(2,1)*t/t3)*r2(i)^3 + mu*D(3,1)*(t^2-t1^2)*t1/t3) / (6*r2(i)^3 + mu*(t^2-t3^2)) - D(1,1) ) ;
        rho(2) = A + mu*B/r2(i)^3 ;
        rho(3) = 1/D0 * ( (6*(D(1,3)*t3/t1 + D(2,3)*t/t1)*r2(i)^3 + mu*D(1,3)*(t^2-t1^3)*t3/t1) / (6*r2(i)^3 + mu*(t^2-t1^2)) - D(3,3) ) ;

        RsatT(1,:) = R(1,:) + rho(1) * Phat(1,:) ;
        RsatT(2,:) = R(2,:) + rho(2) * Phat(2,:) ;
        RsatT(3,:) = R(3,:) + rho(3) * Phat(3,:) ;
        
        v2T = gibbs(RsatT,mu) ; 
        
        [a,e,I,O,W,f] = RV2OE(RsatT(2,:),v2T,mu) ;
        rp = a * (1-e) ;
        
        if ( e >=0 && e < 1) && (rp > Re)
            Rsat = RsatT ;
            V2sat = v2t ;
            OE = [a,e,I,O,W,f] ;
        end
    end
end



function v2 = gibbs(R, mu)
    r(1) = norm(R(1,:)) ;
    r(2) = norm(R(2,:)) ;
    r(3) = norm(R(3,:)) ;
    
    D = cross(R(2,:),R(3,:)) + cross(R(3,:),R(1,:)) + cross(R(1,:),R(2,:)) ;
    N = r(1) * cross(R(2,:),R(3,:)) +  r(2) * cross(R(3,:),R(1,:)) + r(3) * cross(R(1,:),R(2,:)) ;
    S = (r(2)-r(3))*R(1,:) + (r(3)-r(1))*R(2,:) + (r(1)-r(2))*R(3,:) ;
    
    d = norm(D) ;
    n = norm(N) ;

    v2 = 1/r(2) * sqrt( mu / (n*d)) * cross(D,R(2,:)) + sqrt(mu/(n*d)) * S ;
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