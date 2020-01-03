%% AERO 452 Proj1 Kamme
close all
clear all
clc

%% Initial Conditions of the Target A

% Constants
muearth = 398600; %km^3/s^2
rearth = 6378; %km
Tcp = 10*24*60*60; %time to capture in seconds

eccA = 0; % assume circular starting orbit 
incA = 0; % degrees for now
RAANA = deg2rad(115.4157); % degrees for now 
perA = deg2rad(75.9974); % degrees for now
thetaA  = deg2rad(190.4478); % true anomaly
TA = 23.94*60*60; % period in seconds
nA = ((2*pi)/TA); %mean motion of the target
aA = (muearth/(nA^2))^(1/3); %semi major axis
hA = sqrt(aA*muearth*(1-eccA^2));


[rinitialA,vinitialA] = coes2rv(hA,incA,RAANA,eccA,perA,thetaA,muearth);



%% Place the Chaser

% need to place the chaser in the same orbit such that it is 100 km away
% Calculate the difference in true anomaly using chord length

eccB = 0; 
incB = 0;
RAANB = deg2rad(115.4157); % degrees for now 
perB = deg2rad(75.9974); % degrees for now
TB = 23.94*60*60; % period in seconds
nB = ((2*pi)/TB); %mean motion of the target
aB = (muearth/(nB^2))^(1/3); %semi major axis
hB = sqrt(aB*muearth*(1-eccB^2));


% Chord length calculation
% dtheta = asin(100/norm(rinitialA)); % how far ahead in true anomaly the chaser is
% thetaB = thetaA + dtheta; % initial true anomaly of the chaser. Placed ahead of the target


% [rinitialB,vinitialB] = coes2rv(hB,incB,RAANB,eccB,perB,thetaB,muearth);
rHill = [0;100;0]; %km
vHill = [0;0;0];

[rinitialB,vinitialB] =  Hill2ECI_Vectorized(rinitialA,vinitialA,rHill,vHill);


%% Hop Closer--------------------------------------------------------------

% Before Burn
rHill = rHill;
vHill = vHill;

% Burn calc
dVY = (rHill(2) - 40)*(nA/(6*pi));
vHill = vHill + [0;dVY;0];

% After Burn
[rHillF,vHillF] = CWHPropagator(rHill,vHill,nA,TA);

rHillplot(:,1) = rHill;
vHillplot(:,1) = vHill;
rAECI(:,1) = rinitialA;
vAECI(:,1) = vinitialA;
rBECI(:,1) = rinitialB;
vBECI(:,1) = vinitialB;


t = linspace(1,TA);
for i = 2:length(t)
    [rHillplot(:,i),vHillplot(:,i)] = CWHPropagator(rHill,vHill,nA,t(i));
    [rAECI(:,i),vAECI(:,i)] = keplerUniversal(rinitialA,vinitialA,t(i),muearth);
    [rBECI(:,i),vBECI(:,i)] = Hill2ECI_Vectorized(rAECI(:,i),vAECI(:,i),rHillplot(:,i),vHillplot(:,i));
end

originx = 0;
originy = 0;

figure
hold on
title(' Hop to 40 km on Vbar')
plot(rHillplot(2,:),rHillplot(1,:),'k')
plot(originx,originy,'s','Color','k','MarkerSize',10,'MarkerFaceColor','g')
plot(rHillplot(2,1),rHillplot(1,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rHillplot(2,length(rHillplot)),rHillplot(1,length(rHillplot)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
xlabel('Vbar (km)')
ylabel('Rbar (km)')
legend('Hop','Target','Chaser Before Hop','Chaser After Hop')

figure
hold on
title('ECI Hop to 40 km on Vbar')
plot(rAECI(1,:),rAECI(2,:))
plot(rBECI(1,:),rBECI(2,:))
plot(rAECI(1,1),rAECI(2,1),'^','Color','g')
plot(rAECI(1,length(rAECI)),rAECI(2,length(rAECI)),'x','Color','g')
plot(rBECI(1,1),rBECI(2,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rBECI(1,length(rBECI)),rBECI(2,length(rBECI)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
legend('Target Orbit','Chaser Orbit','Target Initial Position','Target Final Position','Chaser Initial Position','Chaser final position')
xlabel('km')
ylabel('km')

%% Time and DV check
clc
timetotal = TA;
dVtotal = 2*dVY;

disp("Time so far: "+timetotal/(60*60*24)+" days")
disp("Dv so far: "+dVtotal*1000+" m/s")






%% Football Motion---------------------------------------------------------

% Before Burn
rHill = rHillF;
vHill = vHillF - [0;dVY;0]; % Need to burn out of the hop

% Burn Calc 
dVX = (norm(rHill)*nA)/2;
vHill = vHill + [dVX;0;0];

% After burn
[rHillF,vHillF] = CWHPropagator(rHill,vHill,nA,TA*1.5);

rHillplotball(:,1) = rHill;
vHillplotball(:,1) = vHill;
rAECI(:,1) = rAECI(:,length(rAECI));
vAECI(:,1) = vAECI(:,length(vAECI));
rBECI(:,1) = rBECI(:,length(rBECI));
vBECI(:,1) = vBECI(:,length(vBECI));


t = linspace(1,TA*1.5);
for i = 2:length(t)
    [rHillplotball(:,i),vHillplotball(:,i)] = CWHPropagator(rHill,vHill,nA,t(i));
    [rAECI(:,i),vAECI(:,i)] = keplerUniversal(rAECI(:,1),vAECI(:,1),t(i),muearth);
    [rBECI(:,i),vBECI(:,i)] = Hill2ECI_Vectorized(rAECI(:,i),vAECI(:,i),rHillplot(:,i),vHillplot(:,i));
end

figure
hold on
plot(rHillplotball(2,:),rHillplotball(1,:),'k')
plot(originx,originy,'s','Color','k','MarkerSize',10,'MarkerFaceColor','g')
% plot(rHillplotball(2,1),rHillplotball(1,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
% plot(rHillplotball(2,lengthball(rHillplotball)),rHillplotball(1,length(rHillplotball)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
xlabel('Vbar (km)')
ylabel('Rbar (km)')
legend('Football','Target')%,'Chaser Before Hop','Chaser After Hop')
title('Football Orbit')


figure
hold on
title('ECI Football')
plot(rAECI(1,:),rAECI(2,:))
plot(rBECI(1,:),rBECI(2,:))
plot(rAECI(1,1),rAECI(2,1),'^','Color','g')
plot(rAECI(1,length(rAECI)),rAECI(2,length(rAECI)),'x','Color','g')
plot(rBECI(1,1),rBECI(2,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rBECI(1,length(rBECI)),rBECI(2,length(rBECI)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
legend('Target Orbit','Chaser Orbit','Target Initial Position','Target Final Position','Chaser Initial Position','Chaser final position')
xlabel('km')
ylabel('km')

%% Time and DV check
clc
timetotal = timetotal+1.5*TA;
dVtotal = dVtotal+2*dVX;

disp("Time so far: "+timetotal/(60*60*24)+" days")
disp("Dv so far: "+dVtotal*1000+" m/s")

%% Hop to Vbar Hold at 1 km------------------------------------------------

% Before Burn
rHill = rHillF;
vHill = vHillF + [dVX;0;0];

% Burn Calc
dr0 = rHill;
drf = [0;-1;0]; %km
dv0minus = vHill;
dvfplus = [0;0;0]; %km/s
[PHIrr,PHIrv,PHIvr,PHIvv] = CWmatrix(nA,TA/2);

dv0plus = inv(PHIrv)*(drf - PHIrr*dr0); 
dvleave = abs(dv0plus - dv0minus);

dvfminus = PHIvr*dr0 + PHIvv*dv0plus;
dvarrive = abs(dvfplus - dvfminus);

% After Burn
rHill = rHill;
vHill = vHill + dv0plus;

[rHillF,vHillF] = CWHPropagator(rHill,vHill,nA,TA*.5);

rHillplot1km(:,1) = rHill;
vHillplot1km(:,1) = vHill;
rAECI(:,1) = rAECI(:,length(rAECI));
vAECI(:,1) = vAECI(:,length(vAECI));
rBECI(:,1) = rBECI(:,length(rBECI));
vBECI(:,1) = vBECI(:,length(vBECI));


t = linspace(1,TA*.5);
for i = 2:length(t)
    [rHillplot1km(:,i),vHillplot1km(:,i)] = CWHPropagator(rHill,vHill,nA,t(i));
    [rAECI(:,i),vAECI(:,i)] = keplerUniversal(rAECI(:,1),vAECI(:,1),t(i),muearth);
    [rBECI(:,i),vBECI(:,i)] = Hill2ECI_Vectorized(rAECI(:,i),vAECI(:,i),rHillplot(:,i),vHillplot(:,i));
end

figure
hold on
plot(rHillplot1km(2,:),rHillplot1km(1,:),'k')
plot(originx,originy,'s','Color','k','MarkerSize',10,'MarkerFaceColor','g')
plot(rHillplot1km(2,1),rHillplot1km(1,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rHillplot1km(2,length(rHillplot1km)),rHillplot1km(1,length(rHillplot1km)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
xlabel('Vbar (km)')
ylabel('Rbar (km)')
legend('2 Impulse','Target','Chaser Before Hop','Chaser After Hop')
title('Hop to Vbar hold at 1km')

figure
hold on
title('ECI Hop to 1 km')
plot(rAECI(1,:),rAECI(2,:))
plot(rBECI(1,:),rBECI(2,:))
plot(rAECI(1,1),rAECI(2,1),'^','Color','g')
plot(rAECI(1,length(rAECI)),rAECI(2,length(rAECI)),'x','Color','g')
plot(rBECI(1,1),rBECI(2,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rBECI(1,length(rBECI)),rBECI(2,length(rBECI)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
legend('Target Orbit','Chaser Orbit','Target Initial Position','Target Final Position','Chaser Initial Position','Chaser final position')
xlabel('km')
ylabel('km')

%% Time and DV check
clc
timetotal = timetotal + .5*TA + 1.0025*TA;
dVtotal = dVtotal + norm(dvleave) + norm(dvarrive);

disp("Time so far: "+timetotal/(60*60*24)+" days")
disp("Dv so far: "+dVtotal*1000+" m/s")


%% Hop to Vbar hold at 300m------------------------------------------------

% Before Burn
rHill = rHillF;
vHill = vHillF - dvarrive;

% Burn Calc
dr0 = rHill;
drf = [0;-.3;0]; %km
dv0minus = vHill;
dvfplus = [0;0;0]; %km/s
[PHIrr,PHIrv,PHIvr,PHIvv] = CWmatrix(nA,TA/2);

dv0plus = inv(PHIrv)*(drf - PHIrr*dr0); 
dvleave = abs(dv0plus - dv0minus);

dvfminus = PHIvr*dr0 + PHIvv*dv0plus;
dvarrive = abs(dvfplus - dvfminus);

% After Burn
rHill = rHill;
vHill = vHill + dv0plus;

[rHillF,vHillF] = CWHPropagator(rHill,vHill,nA,TA*.5);

t = linspace(1,TA*.5);
for i = 1:length(t)
    [rHillplot300m(:,i),vHillplot300m(:,i)] = CWHPropagator(rHill,vHill,nA,t(i));
end

figure
hold on
plot(rHillplot300m(2,:),rHillplot300m(1,:),'k')
plot(originx,originy,'s','Color','k','MarkerSize',10,'MarkerFaceColor','g')
plot(rHillplot300m(2,1),rHillplot300m(1,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rHillplot300m(2,length(rHillplot300m)),rHillplot300m(1,length(rHillplot300m)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
xlabel('Vbar (km)')
ylabel('Rbar (km)')
legend('2 Impulse','Target','Chaser Before Hop','Chaser After Hop')
title('Hop to Vbar hold at 300 m')


%% Time and DV check
clc
timetotal = timetotal + .5*TA + 1.0025*TA;
dVtotal = dVtotal + norm(dvleave) + norm(dvarrive);

disp("Time so far: "+timetotal/(60*60*24)+" days")
disp("Dv so far: "+dVtotal*1000+" m/s")

%% Hop to 20m on Rbar------------------------------------------------------

% Before Burn
rHill = rHillF;
vHill = vHillF - dvarrive;

% Burn Calc
dr0 = rHill;
drf = [0.02;0;0]; %km
dv0minus = vHill;
dvfplus = [0;0;0]; %km/s
[PHIrr,PHIrv,PHIvr,PHIvv] = CWmatrix(nA,.05*TA);

dv0plus = inv(PHIrv)*(drf - PHIrr*dr0); 
dvleave = abs(dv0plus - dv0minus);

dvfminus = PHIvr*dr0 + PHIvv*dv0plus;
dvarrive = abs(dvfplus - dvfminus);

% After Burn
rHill = rHill;
vHill = vHill + dv0plus;

[rHillF,vHillF] = CWHPropagator(rHill,vHill,nA,TA*.05);


t = linspace(1,TA*.05);
for i = 1:length(t)
    [rHillplot20m(:,i),vHillplot20m(:,i)] = CWHPropagator(rHill,vHill,nA,t(i));
end


figure
hold on
plot(rHillplot20m(2,:),rHillplot20m(1,:),'k')
plot(originx,originy,'s','Color','k','MarkerSize',10,'MarkerFaceColor','g')
plot(rHillplot20m(2,1),rHillplot20m(1,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rHillplot20m(2,length(rHillplot20m)),rHillplot20m(1,length(rHillplot20m)),'s','Color','k','MarkerSize',10,'MarkerFaceColor','m')
xlabel('Vbar (km)')
ylabel('Rbar (km)')
legend('2 Impulse','Target','Chaser Before Hop','Chaser After Hop')
title('Hop to Rbar at 20 m')

% DV for vbar hold
time20 = 24*60*60;
a = [-3*nA^2*.02;0;0];
dv20 = a(1)*time20;

%% Time and DV check
clc
timetotal = timetotal + time20 + .05*TA;
dVtotal = dVtotal + norm(dvleave) + norm(dvarrive)+dv20;

disp("Time so far: "+timetotal/(60*60*24)+" days")
disp("Dv so far: "+dVtotal*1000+" m/s")





%% Hop to Rendezvous

% Before Burn
rHill = rHillF;
vHill = dvfplus;

% Try to optimize this hop
timevec = linspace(1,TA,TA);

for i = 1:length(timevec)
% Burn Calc
dr0 = rHill;
drf = [0;0;0]; %km
dv0minus = vHill;
dvfplus = [0;0;0]; %km/s
[PHIrr,PHIrv,PHIvr,PHIvv] = CWmatrix(nA,timevec(i));

dv0plus(:,i) = inv(PHIrv)*(drf - PHIrr*dr0); 
dvleave = abs(dv0plus - dv0minus);

dvfminus = PHIvr*dr0 + PHIvv*dv0plus;
dvarrive = abs(dvfplus - dvfminus);

dv0plusAccumulator(i) = norm(dv0plus(:,i));
end



[M,I] = min(dv0plusAccumulator);

% After Burn
rHill = rHill;
vHill = vHill + dv0plus(:,I);

[rHillF,vHillF] = CWHPropagator(rHill,vHill,nA,I);

t = linspace(1,I);
for i = 1:length(t)
    [rHillplotfinal(:,i),vHillplotfinal(:,i)] = CWHPropagator(rHill,vHill,nA,t(i));
end

figure
hold on
plot(rHillplotfinal(2,:),rHillplotfinal(1,:),'k')
plot(originx,originy,'s','Color','k','MarkerSize',10,'MarkerFaceColor','g')
plot(rHillplotfinal(2,1),rHillplotfinal(1,1),'s','Color','k','MarkerSize',10,'MarkerFaceColor','r')
plot(rHillplotfinal(2,length(rHillplotfinal)),rHillplotfinal(1,length(rHillplotfinal)),'^','Color','k','MarkerSize',5)
xlabel('Vbar (km)')
ylabel('Rbar (km)')
legend('2 Impulse','Target','Chaser Before Hop','Chaser After Hop')
title('Hop to Rendezvous')


%% Time and DV check
clc
timetotal = timetotal + timevec(I) + 60*60*24;
dVtotal = dVtotal+2*norm(dv0plus(:,I));

disp("Time so far: "+timetotal/(60*60*24)+" days")
disp("Dv so far: "+dVtotal*1000+" m/s")


%% Functions

% Coes to R and V
function [R,V] = coes2rv(h,inc,RAAN,e,per,theta,muearth)
% h [km^2/s] Specific angular momentum
% i [rad] Inclination
% RAAN [rad] Right ascension (RA) of the ascending node
% e Eccentricity
% per [rad] Argument of perigee
% theta [rad] True anomaly
% muearth = 398600; Earth’s gravitational parameter [km^3/s^2]

% State Vectors in Perifocal coordinates
rx = h^2/muearth*(1/(1 + e*cos(theta)))*[cos(theta);sin(theta);0];
vx = muearth/h*[-sin(theta); (e +cos(theta));0];

% Direction cosine matrix
DCM = [cos(per), sin(per),0;-sin(per),cos(per),0;0,0,1]*...
 [1,0,0;0,cos(inc),sin(inc);0,-sin(inc),cos(inc)]*...
 [cos(RAAN), sin(RAAN),0;-sin(RAAN),cos(RAAN),0;0,0,1];

% Transformation Matrix
Dcm = inv(DCM);

% ECI R
R = Dcm*rx;

% ECI V
V = Dcm*vx;

end

% Vallado Universial Variable Functions
function [r,v] = keplerUniversal(r0,v0,t,mu)
%input vectors as COLUMNS

%Purpose:
%Most effecient way to propagate any type of two body orbit using kepler's
%equations.
%-------------------------------------------------------------------------%
%                                                                         %
% Inputs:                                                                 %
%--------                                                                  
%r_ECI                  [3 x N]                         Position Vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%v_ECI                  [3 x N]                         Velocity vector in
%                                                       ECI coordinate
%                                                       frame of reference
%
%t                      [1 x N]                         time vector in
%                                                       seconds
%
%mu                     double                          Gravitational Constant
%                                                       Defaults to Earth if
%                                                       not specified
% Outputs:
%---------                                                                %
%r_ECI                  [3 x N]                         Final position
%                                                       vector in ECI
%
%v_ECI                  [3 x N]                         Final velocity
%                                                       vector in ECI
%--------------------------------------------------------------------------
% Programmed by Darin Koblick 03-04-2012                                  %
%-------------------------------------------------------------------------- 
if ~exist('mu','var'); mu = 398600.4418; end
tol = 1e-9;
v0Mag = sqrt(sum(v0.^2,1));  r0Mag = sqrt(sum(r0.^2,1));
alpha = -(v0Mag.^2)./mu + 2./r0Mag; 
% Compute initial guess (X0) for Newton's Method
X0 = NaN(size(t));
%Check if there are any Eliptic/Circular orbits
idx = alpha > 0.000001;
if any(idx)
    X0(idx) = sqrt(mu).*t(idx).*alpha(idx); 
end
%Check if there are any Parabolic orbits
idx = abs(alpha) < 0.000001;
if any(idx)
   h = cross(r0(:,idx),v0(:,idx)); hMag = sqrt(sum(h.^2,1));
   p = (hMag.^2)./mu; s = acot(3.*sqrt(mu./(p.^3)).*t(idx))./2;
   w = atan(tan(s).^(1/3)); X0(idx) = sqrt(p).*2.*cot(2.*w);
end
%Check if there are any Hyperbolic orbits
idx = alpha < -0.000001;
if any(idx)
   a = 1./alpha(idx);
   X0(idx) = sign(t(idx)).*sqrt(-a).*...
       log(-2.*mu.*alpha(idx).*t(idx)./ ...
       (dot(r0(:,idx),v0(:,idx))+sign(t(idx)).*sqrt(-mu.*a).*...
       (1-r0Mag(idx).*alpha(idx))));
end
% Newton's Method to converge on solution
% Declare Constants that do not need to be computed within the while loop
err = Inf;
dr0v0Smu = dot(r0,v0)./sqrt(mu);
Smut = sqrt(mu).*t;
while any(abs(err) > tol)
    X02 = X0.^2;
    X03 = X02.*X0;
    psi = X02.*alpha;
    [c2,c3] = c2c3(psi);
    X0tOmPsiC3 = X0.*(1-psi.*c3);
    X02tC2 = X02.*c2;
    r = X02tC2 + dr0v0Smu.*X0tOmPsiC3 + r0Mag.*(1-psi.*c2);
    Xn = X0 + (Smut-X03.*c3-dr0v0Smu.*X02tC2-r0Mag.*X0tOmPsiC3)./r;
    err = Xn-X0; X0 = Xn;
end
f = 1 - (Xn.^2).*c2./r0Mag; g = t - (Xn.^3).*c3./sqrt(mu);
gdot = 1 - c2.*(Xn.^2)./r; fdot = Xn.*(psi.*c3-1).*sqrt(mu)./(r.*r0Mag);
r = bsxfun(@times,f,r0) + bsxfun(@times,g,v0);
v = bsxfun(@times,fdot,r0) + bsxfun(@times,gdot,v0);
% Ensure Solution Integrity
%idx = round((f.*gdot - fdot.*g)./tol).*tol ~= 1; r(:,idx) = NaN; v(:,idx) = NaN;
end

function [c2,c3] = c2c3(psi)
%Vallado pg. 71 Algorithm 1
c2 = NaN(size(psi));
c3 = NaN(size(psi));
idx = psi > 1e-6;
if any(idx)
    c2(idx) = (1-cos(sqrt(psi(idx))))./psi(idx);
    c3(idx) = (sqrt(psi(idx))-sin(sqrt(psi(idx))))./sqrt(psi(idx).^3);
end
idx = psi < -1e-6;
if any(idx)
    c2(idx) = (1 - cosh(sqrt(-psi(idx))))./psi(idx);
    c3(idx) = (sinh(sqrt(-psi(idx)))-sqrt(-psi(idx)))./sqrt(-psi(idx).^3);
end
idx = abs(psi) <= 1e-6;
if any(idx)
    c2(idx) = 0.5;
    c3(idx) = 1/6;
end
end

% CW matrices
function [PHIrr,PHIrv,PHIvr,PHIvv] = CWmatrix(n,t)
PHIrr = [(4-3*cos(n*t)) 0 0;
         6*(sin(n*t)-n*t) 1 0;
         0 0 cos(n*t)];
         
PHIrv = [(1/n)*sin(n*t) (2/n)*(1-cos(n*t)) 0;
         (2/n)*(cos(n*t)-1) (1/n)*(4*sin(n*t)-3*n*t) 0;
         0 0 (1/n)*sin(n*t)];
         
PHIvr = [3*n*sin(n*t) 0 0;
         6*n*(cos(n*t)-1) 0 0;
         0 0 (-n*sin(n*t))];
     
PHIvv = [cos(n*t) 2*sin(n*t) 0;
         -2*sin(n*t) (4*cos(n*t)-3) 0;
         0 0 cos(n*t)];
end

% Calculation of relative things using actual relative motion equations
function [rrelLVLH,vrelLVLH,arelLVLH] = relthings(rAvec,vAvec,rBvec,vBvec)
% Vector inputs must be COLUMNS in ECI
% A is target and B is the chaser
muearth = 398600;


hAvec = cross(rAvec,vAvec);

ihatA = ((rAvec)/norm(rAvec));

khatA = (hAvec)/norm(hAvec);

jhatA = cross(khatA,ihatA);

Qxx = [ihatA';jhatA';khatA']; %orthogonal transformation matrix

% Calculation of relative position in ECI
rrel = rBvec - rAvec;


angvelA = ((hAvec)/(norm(rAvec))^2); %angular velocity of the comoving A frame

angaccA = (-2)*((dot(vAvec,rAvec))/(norm(rAvec))^2)*angvelA; %angular acceleration of the comoving target A frame

% Calculation of absolute accelerations
aA = (-muearth)*(rAvec/(norm(rAvec))^3); % km/s^2 absolute acceleration of S/C A
aB = (-muearth)*(rBvec/(norm(rBvec))^3); % km/s^2 absolute accleration of S/C B



% Calculation of relative velocity
vrel = vBvec - vAvec - cross(angvelA,rrel);

% Calculation of relative accleration
arel = aB - aA - cross(angaccA,rrel) - cross(angvelA,cross(angvelA,rrel)) - 2*cross(angvelA,vrel);

% Rotating relative position/velocity/acceleration into comoving frame of target S/C A
rrelLVLH = (Qxx*rrel);
vrelLVLH = (Qxx*vrel);
arelLVLH = (Qxx*arel);
end

% Vallado Hill Stuff
function [rHill,vHill] = ECI2Hill_Vectorized(rTgt,vTgt,rChase,vChase)
% Purpose:
% Convert those position (ECI) and velocity (ECI) into Hill's reference
% frame using both the target and the chaser position/velocity data
%
% Inputs:
%rTgt                       [3 x N]                 ECI Position vector of
%                                                   reference frame (km)
%
%vTgt                       [3 x N]                 ECI Velocity vector of
%                                                   reference frame (km/s)
%rChase                     [3 x N]
%
%vChase                     [3 x N]
%
% Outputs:
%rHill                      [3 x N]                 Hill's relative
%                                                   position vector (km)
%
%vHill                      [3 x N]                 Hill's relative
%                                                   velocity vector (km/s)
% References:
% Vallado 2007.
% Programed by Darin C Koblick 11/30/2012
% Begin Code Sequence
%Declare Local Functions
matrixMultiply = @(x,y)permute(cat(2,sum(permute(x(1,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(2,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(3,:,:),[3 2 1]).*permute(y,[2 1]),2)),[2 1]);
rTgtMag = sqrt(sum(rTgt.^2,1));
rChaseMag = sqrt(sum(rChase.^2,1));
vTgtMag = sqrt(sum(vTgt.^2,1));
%Determine the RSW transformation matrix from the target frame of reference
RSW = ECI2RSW(rTgt,vTgt);
%Use RSW rotation matrix to convert rChase and vChase to RSW
r_Chase_RSW = matrixMultiply(RSW,rChase);
v_Chase_RSW = matrixMultiply(RSW,vChase);
%Find Rotation angles to go from target to interceptor
phi_chase = asin(r_Chase_RSW(3,:)./rChaseMag);
lambda_chase = atan2(r_Chase_RSW(2,:),r_Chase_RSW(1,:));
CPC = cos(phi_chase);     SPC = sin(phi_chase);
SLC = sin(lambda_chase);  CLC = cos(lambda_chase);
%Find Position component rotations
rHill = cat(1,rChaseMag-rTgtMag, ...
              lambda_chase.*rTgtMag, ...
              phi_chase.*rTgtMag);
%Find the rotation matrix RSW->SEZ of chaser
RSW_SEZ = zeros(3,3,size(rTgtMag,2));
RSW_SEZ(1,1,:) = SPC.*CLC;  RSW_SEZ(1,2,:) = SPC.*SLC;  RSW_SEZ(1,3,:) = -CPC;
RSW_SEZ(2,1,:) = -SLC;  RSW_SEZ(2,2,:) = CLC;
RSW_SEZ(3,1,:) = CPC.*CLC;  RSW_SEZ(3,2,:) = CPC.*SLC;  RSW_SEZ(3,3,:) = SPC;
%Find the velocity component of positions using the angular rates in SEZ frame
v_Chase_SEZ = matrixMultiply(RSW_SEZ,v_Chase_RSW);
vHill = cat(1,v_Chase_SEZ(3,:), ...
              rTgtMag.*(v_Chase_SEZ(2,:)./(rChaseMag.*CPC)-vTgtMag./rTgtMag), ...
              -rTgtMag.*v_Chase_SEZ(1,:)./rChaseMag);
end

function [T,rRSW,vRSW] = ECI2RSW(rECI,vECI)
% Purpose:
%Convert ECI Coordinates to RSW Coordinates, also, return the
%transformation matrix T in which to take a given set of coordinates in ECI
%and convert them using the same RSW reference frame.
%
% Inputs:
% rECI              [3 x N]                     ECI position Coordinates in
%                                               km
%
% vECI              [3 x N]                     ECI velocity Coordinates in
%                                               km/s
%
% Outputs:
% T                 [3 x 3 x N]                 Transformation matrix
%                                               necessary to go from
%                                               rECI -> rRSW
%
% rRSW              [3 x N]                     RSW Position Coordinates
%                                               km
%
% vRSW              [3 x N]                     RSW Velocity Coordinates
%                                               km/s
%
% References:
% Vallado pg. 173
% Vallado rv2rsw.m code dated 06/09/2002
%
%Programmed by: Darin C Koblick                 11/29/2012
% Begin Code Sequence
if nargin == 0
    rECI =  repmat([6968.1363,1,2]',[1 10]);
    vECI =  repmat([3,7.90536615282099,4]',[1 10]);
    [T,rRSW,vRSW] = ECI2RSW(rECI,vECI);
    return;
end
%Declared Internal function set:
unitv = @(x)bsxfun(@rdivide,x,sqrt(sum(x.^2,1)));
matrixMultiply = @(x,y)permute(cat(2,sum(permute(x(1,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(2,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(3,:,:),[3 2 1]).*permute(y,[2 1]),2)),[2 1]);
%Find the Radial component of the RIC position vector
rvec = unitv(rECI);
%Find the cross-track component of the RIC position vector
wvec = unitv(cross(rECI,vECI));
%Find the along-track component of the RIC position vector
svec = unitv(cross(wvec,rvec));
%Create the transformation matrix from ECI to RSW
T = NaN(3,3,size(rECI,2));
T(1,1,:) = rvec(1,:); T(1,2,:) = rvec(2,:); T(1,3,:) = rvec(3,:);
T(2,1,:) = svec(1,:); T(2,2,:) = svec(2,:); T(2,3,:) = svec(3,:);
T(3,1,:) = wvec(1,:); T(3,2,:) = wvec(2,:); T(3,3,:) = wvec(3,:);
%Find the position and velocity vectors in the RSW reference frame!
rRSW = matrixMultiply(T,rECI);
vRSW = matrixMultiply(T,vECI);
end

function [rHill,vHill] = CWHPropagator(rHillInit,vHillInit,omega,t)
% Purpose:
% Take initial position and velocity coordinates in the Hill reference frame
% and propagate them using the Clohessy-Wiltshire Hill Linearize equation
% of motion.
%
% Inputs:
%rHillInit                  [3 x 1]                 Hill Position vector
%                                                   (km) / (m)
%
%vHillInit                  [3 x 1]                 Hill Velocity vector of
%                                                   (km/s) / (m/s)
%
%omega                       double                 Orbital Angular Rate
%                                                   of the target
%                                                   (rad/s)
%                                                   Should be close to
%                                                   circular for linear propagation
%                                                   error to be low.
%
%t                          [1 x N]                 Propagation Time in
%                                                   seconds
%                                                   
%
%
%
% Outputs:
%rHill                       [3 x N]                Propagated Hill
%                                                   Position vector (km) /
%                                                   (m/s)
%
%vHill                       [3 x N]                Propagated Hill
%                                                   Velocity vector (km/s)
%                                                   / (m/s)
%
%
% References:
%
% Programed by Darin C Koblick 11/30/2012
% Begin Code Sequence
x0 = rHillInit(1,:); y0 = rHillInit(2,:); z0 = rHillInit(3,:);
x0dot = vHillInit(1,:); y0dot = vHillInit(2,:); z0dot = vHillInit(3,:);
rHill = [(x0dot./omega).*sin(omega.*t)-(3.*x0+2.*y0dot./omega).*cos(omega.*t)+(4.*x0+2.*y0dot./omega)
        (6.*x0+4.*y0dot./omega).*sin(omega.*t)+2.*(x0dot./omega).*cos(omega.*t)-(6.*omega.*x0+3.*y0dot).*t+(y0-2.*x0dot./omega)
        z0.*cos(omega.*t)+(z0dot./omega).*sin(omega.*t)];
vHill = [x0dot.*cos(omega.*t)+(3.*omega.*x0+2.*y0dot).*sin(omega.*t)
        (6.*omega.*x0 + 4.*y0dot).*cos(omega.*t) - 2.*x0dot.*sin(omega.*t)-(6.*omega.*x0 + 3.*y0dot)
        -z0.*omega.*sin(omega.*t)+z0dot.*cos(omega.*t)];
end

function [rInt,vInt] = Hill2ECI_Vectorized(rTgt,vTgt,rHill,vHill)
% Purpose:
% Convert those position (rHill) and velocity (vHill) values back into an
% ECI coordinate frame of reference using the reference satellite
% (rTgt,vTgt) position and velocity data.
%
% Inputs:
%rTgt                       [3 x N]                 ECI Position vector of
%                                                   reference frame (km)
%
%vTgt                       [3 x N]                 ECI Velocity vector of
%                                                   reference frame (km/s)
%
%rHill                      [3 x N]                 Hill's relative
%                                                   position vector (km)
%
%vHill                      [3 x N]                 Hill's relative
%                                                   velocity vector (km/s)
%
%
%
% Outputs:
%rInt                       [3 x N]
%
%vInt                       [3 x N]
%
%
% References:
% Vallado 2007.
% Programed by Darin C Koblick 11/30/2012
% Begin Code Sequence
%Declare Local Functions
rTgtMag = sqrt(sum(rTgt.^2,1));
vTgtMag = sqrt(sum(vTgt.^2,1));
matrixMultiply = @(x,y)permute(cat(2,sum(permute(x(1,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(2,:,:),[3 2 1]).*permute(y,[2 1]),2), ...
                                     sum(permute(x(3,:,:),[3 2 1]).*permute(y,[2 1]),2)),[2 1]);
%Find the RSW matrix from the target ECI positions
RSW = ECI2RSW(rTgt,vTgt); rIntMag = rTgtMag + rHill(1,:);
%Compute rotation angles to go from tgt to interceptor
lambda_int = rHill(2,:)./rTgtMag; phi_int = sin(rHill(3,:)./rTgtMag);
CLI = cos(lambda_int); SLI = sin(lambda_int); CPI = cos(phi_int); SPI = sin(phi_int);
%find rotation matrix to go from rsw to SEZ of inerceptor
RSW_SEZ = zeros(3,3,size(rTgt,2));
RSW_SEZ(1,1,:) = SPI.*CLI;  RSW_SEZ(1,2,:) = SPI.*SLI; RSW_SEZ(1,3,:) = -CPI;
RSW_SEZ(2,1,:) = -SLI;      RSW_SEZ(2,2,:) = CLI;      RSW_SEZ(3,1,:) = CPI.*CLI;
                            RSW_SEZ(3,2,:) = CPI.*SLI; RSW_SEZ(3,3,:) = SPI;
%Find velocity component positions by using angular rates in SEZ frame
vIntSEZ = cat(1,-rIntMag.*vHill(3,:)./rTgtMag, ...
                 rIntMag.*(vHill(2,:)./rTgtMag + vTgtMag./rTgtMag).*CPI, ...
                 vHill(1,:));
vInt = matrixMultiply(permute(RSW,[2 1 3]), ...
       matrixMultiply(permute(RSW_SEZ,[2 1 3]), ...
       vIntSEZ));
%Find the position components
rIntRSW = bsxfun(@times,rIntMag,cat(1,CPI.*CLI, ...
                                      CPI.*SLI, ...
                                      SPI));
rInt = matrixMultiply(permute(RSW,[2 1 3]),rIntRSW);      
end