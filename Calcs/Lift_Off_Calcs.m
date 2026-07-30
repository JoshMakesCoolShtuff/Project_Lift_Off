W = 8.456;
T = 1*4.44822;
CL = 0.4; %@ 0 AOA
S = 0.079;



[~,~,~,rho] = atmosisa(0);


%min takeoff speed

Vmin = sqrt(W/(0.5*rho*S*CL))

