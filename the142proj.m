clearvars
close all
%notes in ECI system, for 30 lat and long area! So that where you need to

R = 6371; % earth radius in km
latspacing = 10; 
lonspacing = 20; 
% lines of longitude: 
[lon1,lat1] = meshgrid(-180:lonspacing:180,linspace(-90,90,300)); 
[x1,y1,z1] = sph2cart(lon1*pi/180,lat1*pi/180,R); 
plot3(x1,y1,z1,'-','color',0.5*[1 1 1])
hold on
[X,Y,Z] = sphere(100);
surf(X*R*.99,Y*R*.99,Z*R*.99)
% lines of latitude: 
[lat2,lon2] = meshgrid(-90:latspacing:90,linspace(-180,180,300)); 
[x2,y2,z2] = sph2cart(lon2*pi/180,lat2*pi/180,R); 
plot3(x2,y2,z2,'-','color',0.5*[1 1 1])

xlim([0,R+1000])
ylim([0,R+1000])
zlim([0,R+1000])

%30 lat and long
[x3,y3,z3] = sph2cart([0:1:90].*pi./180,30.*pi./180.*ones(1,91),R);
plot3(x3,y3,z3,'-','color','red',LineWidth=3)
[x4,y4,z4] = sph2cart(30.*pi./180.*ones(1,91),[0:1:90].*pi./180,R);
plot3(x4,y4,z4,'-','color','red',LineWidth=3)
%0 lat and long
[x5,y5,z5] = sph2cart([0:1:90].*pi./180,0.*pi./180.*ones(1,91),R);
plot3(x5,y5,z5,'-','color','blue',LineWidth=3)
[x6,y6,z6] = sph2cart(0.*pi./180.*ones(1,91),[0:1:90].*pi./180,R);
plot3(x6,y6,z6,'-','color','blue',LineWidth=3)


%zoom in on it
% xlim([0,8600])
% ylim([0,3500])
% zlim([0,3400])



%% Choose satelite locations
P1 = [7000 1300 2800];
P2 = [7200 500 1100];
P3 = [7300 2500 1500];
scatter3(P1(1),P1(2),P1(3),'o','blue','LineWidth',4)
scatter3(P2(1),P2(2),P2(3),'o','green','LineWidth',4)
scatter3(P3(1),P3(2),P3(3),'o','magenta','LineWidth',4)

%initial state vector
[x0,y0,z0] = sph2cart(28.*pi./180,5.*pi./180,R+100);
x0 = [x0 y0 z0]';
v0 = zeros(0,3);

%dynamics model, constant linear velocty
v = [5; -30; 25];
v0 = v;
x = zeros(3,31);
x(:,1) = x0;
%determine simulation time
dt = 1;
simTime = 30; % seconds
t = 0:dt:simTime;
%t = linspace(0,100,50);
delT = dt; %time step
%propogate dynamics through time
for i = 2:length(t)
    x(:,i) = x(:,i-1)+v*delT;
end

%% Calculate pseudoranges
%Using dynamics, get what the new postition should be with some added
%error.
%Find the pseudoranges via the simulation location.

plot3(x(1,:),x(2,:),x(3,:),'cyan',"LineWidth",3);

% 
% for i = 1:length(t)
%     pr1(i) = norm([P1(1)-x(1,i),P1(2)-x(2,i),P1(3)-x(3,i)]);
%     pr2(i) = norm([P2(1)-x(1,i),P2(2)-x(2,i),P2(3)-x(3,i)]);
%     pr3(i) = norm([P3(1)-x(1,i),P3(2)-x(2,i),P3(3)-x(3,i)]);
% 
% end
% 
% %% Estimate Receiver Position from Pseudoranges and Satellite Positions
% %Trilateration
% 
% Cord1 = ones(3,length(t));
% Cord2 = ones(3,length(t));
% for i = 1:length(t)
%     [Cord1(:,i),Cord2(:,i)] = trilat(P1,P2,P3,pr1(i),pr2(i),pr3(i));
% end
% 
% plot3(Cord1(:,1),Cord1(:,2),Cord1(:,3));
% plot3(Cord2(:,1),Cord2(:,2),Cord2(:,3));
% legend("actual","cord1","cord2")
% disp("end")
% 
% function [K1,K2] = trilat(P1, P2, P3, r1, r2, r3)
% 
%   p1 = [0, 0, 0];
%   p2 = [P2(1) - P1(1), P2(2) - P1(2), P2(3) - P1(3)];
%   p3 = [P3(1) - P1(1), P3(2) - P1(2), P3(3) - P1(3)];
%   v1 = p2 - p1;
%   v2 = p3 - p1;
% 
%   Xn = (v1)/norm(v1);
% 
%   tmp = cross(v1, v2);
% 
%   Zn = tmp/norm(tmp);
% 
%   Yn = cross(Xn, Zn);
% 
%   i = dot(Xn, v2);
%   d = dot(Xn, v1);
%   j = dot(Yn, v2);
% 
%   X = ((r1^2)-(r2^2)+(d^2))/(2*d);
%   Y = (((r1^2)-(r3^2)+(i^2)+(j^2))/(2*j))-((i/j)*(X));
%   Z1 = sqrt(max(0, r1^2-X^2-Y^2));
%   Z2 = -Z1;
% 
%   K1 = P1 + X * Xn + Y * Yn + Z1 * Zn;
%   K2 = P1 + X * Xn + Y * Yn + Z2 * Zn;
% end

%% Generic Kalman filter

%% Parameters
n = 6;      % number of dimensions
dt = 1;     % seconds
simTime = 30; % seconds
tspan = 0:dt:simTime;
numsteps = length(tspan);
%R = 10000; % radius of earth

% initial conditions
r0 = [x0(1), x0(2), x0(3), 5, -30, 25];

% motion model dynamics
A = [1 dt 0 0  0 0;
     0 1  0 0  0 0;
     0 0  1 dt 0 0;
     0 0  0 1  0 0;
     0 0  0 0  1 dt; 
     0 0  0 0  0 1];

% measurement dynamics
% H = [1 0 0 0 0 0;
%      0 0 1 0 0 0;
%      0 0 0 0 1 0;];
H = [1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0];

% process noise covariance
Cw = 0.05;
B = diag(ones(1,n));
% measurement noise covariance
Cv = 0.01;


% initialize vectors
y = NaN(3, numsteps);
M = zeros(n,n,numsteps);
P = NaN(n,n,numsteps);
r_true = NaN(n,numsteps);
r_estimate = NaN(n,numsteps);
r_nominal = NaN(n,numsteps);

% store intial conditions
r_nominal(:,1) = r0;
r_true(:,1) = r0;
r_estimate(:,1) = r0;
M(:,:,1) = B*Cw*B';


for i = 1:numsteps-1 % Kalman Loop

    %%% DYNAMICS UPDATE %%%
    % get latest apriori
    xk_apriori = A*r_estimate(:,i);

    %%% OBSERVATION UPDATE %%%
    % calculate satellite observation y
    yk = H*r_true(:,i) + randn(3,1)*sqrt(Cv); %r_true = z_k?
    y(:,i) = yk;
    
    % Kalman gain
    Mk = M(:,:,i); 
    Pk = inv(inv(Mk) + H'*inv(Cv)*H);
    Kalman_gain = Pk*H'*inv(Cv);
    
    % find aposteriori
    xk_apost = xk_apriori + Kalman_gain*(yk-H*xk_apriori);
    
    %%% ITERATION UPDATES %%% 

    % update states
    r_nominal(:,i+1) = A*r_nominal(:,i);
    r_true(:,i+1) = A*r_true(:,i) + B*randn(n,1)*sqrt(Cw);
    r_estimate(:,i+1) = xk_apost;

    % update covariances
    M(:,:,i+1) = A*Pk*A' + B*Cw*B';

end

plot3(r_true(1,:),r_true(2,:),r_true(3,:),'red',"LineWidth",3)



disp("done")