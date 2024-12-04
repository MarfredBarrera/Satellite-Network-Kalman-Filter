%% Generic Kalman filter

%% Parameters
n = 6;      % number of dimensions
dt = 1;     % seconds
simTime = 30; % seconds
tspan = 0:dt:simTime;
numsteps = length(tspan);
R = 10000; % radius of earth

% initial conditions
r0 = [0, 1, 0, 0, 0, 0];

% motion model dynamics
A = [1 dt 0 0  0 0;
     0 1  0 0  0 0;
     0 0  1 dt 0 0;
     0 0  0 1  0 0;
     0 0  0 0  1 dt; 
     0 0  0 0  0 1];

% measurement dynamics
H = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0;];

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
    yk = H*r_true(:,i) + randn(3,1)*sqrt(Cv);
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






    











