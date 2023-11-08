clc
clear
close all

load('IODMeasurements.mat')
load('IODMeasurements2.mat')

l = length(AZIMUTH) ;
LOS = zeros(l,3) ;

for i = 1:l
    LOS(i,1) = cosd(ELEVATION(i)) * cosd(AZIMUTH(i)) ;
    LOS(i,2) = cosd(ELEVATION(i)) * sind(AZIMUTH(i)) ;
    LOS(i,3) = sind(ELEVATION(i)) ;
end

RSAT = zeros(3,3,l/3) ;
for i = 1:l/3
    RSAT(:,:,i) = gauss(TIMES((i*3-2:i*3)) ,RSITES((i*3-2:i*3),:) ,LOS((i*3-2:i*3),:),MU,Re)
    
end 



function Rsat = gauss(time,R,Phat,mu,Re)
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
    r2 = r2(r2==real(r2)) ;
    %rejects values that are not greater than the radius of the earth
    r2 = r2(r2>Re) 
    
    rho(1) = 1/D0 * ( (6*(D(3,1)*t1/t3 + D(2,1)*t/t3) + mu*D(3,1)*(t^2-t1^2)*t1/t3) / (6*r2^3 + mu*(t^2-t3^2)) - D(1,1) ) ;
    rho(2) = A + mu*B/r2^3 ;
    rho(3) = 1/D0 * ( (6*(D(1,3)*t3/t1 + D(2,3)*t/t1) + mu*D(1,3)*(t^2-t1^3)*t3/t1) / (6*r2^3 + mu*(t^2-t1^2)) - D(3,3) ) ;

    Rsat(1,:) = R(1,:) + rho(1) * p(1,:) ;
    Rsat(2,:) = R(2,:) + rho(2) * p(2,:) ;
    Rsat(3,:) = R(3,:) + rho(3) * p(3,:) 
end



function v2 = gibbs(r, mu)
    
    

end
