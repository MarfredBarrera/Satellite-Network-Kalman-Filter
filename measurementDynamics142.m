%% Measurement dynamics 

%% Given
% Position of each satellite

P1 = [7000 1300 2800]';
P2 = [7200 500 1100]';
P3 = [7300 2500 1500]';
R_sat = [P1 P2 P3];

% Assuming we know the position of the aircraft the moment that its GPS
% goes offline, we have initial position of aircraft.
R_ac = [6800; 1000; 1000];

% speed of light
c = 299792; %km/s

% time emitted
t = 0;

H = zeros(3,3);
for i = 1:3

    d(:,i) = abs(R_sat(:,i) - R_ac); % [3x3] matrix which contains distance decomposition [d1(x;y;z) d2(x;y;z) d3(x;y;z)] between each satellite and the aircraft

    % denom(i) = sqrt( (R_ac(1)-R_sat(1,i))^2 + (R_ac(2)-R_sat(2,i))^2 + (R_ac(3)-R_sat(3,i))^2 ) ; 
    denom(i) = sqrt( (d(1,i)^2 + d(2,i)^2 + d(3,i)^2)); % [1x3] matrix of magnitude of distance from each of the satellites

    % H(i,:) = [ (R_sat(1,i)-R_ac(1))/denom(i) (R_sat(2,i)-R_ac(2))/denom(i)  (R_sat(3,i)-R_ac(3))/denom(i) ];
    H(i,:) = [ d(1,i)/denom(i) d(2,i)/denom(i)  d(3,i)/denom(i) ];



end


% H = H./c;
% 
% H(:,4) = 1





