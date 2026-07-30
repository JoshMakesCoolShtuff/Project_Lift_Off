function xdot = CessnaEOM(x,u)
    %% States and Inputs
    V     = x(1);
    alpha = x(2);
    q     = x(3);
    theta = x(4);
    p     = x(5);
    h     = x(6);
    
    el = u(1);
    th = u(2);
    
    %% Cessna Constants
    
    g   = 32.17;
    rho = 0.002377;
    
    W    = 2650;
    S    = 174;
    cbar = 4.9;
    Jz   = 1346;
    m    = W/g;
    
    CL0   = 0.307;
    CLa   = 4.41;
    CL_el = 0.43;
    
    CDM  = 0.0223;
    k    = 0.0554;
    CLDM = 0;
    
    CM_R0  = 0.04;
    CM_Ra  = -0.613;
    CM_Rel = -1.122;
    
    CM_adot = -7.27;
    CM_q    = -12.4;
    CL_adot = 1.7;
    CL_q    = 3.9;
    
    eta = 0.7;
    
    epsilon = 0;
    eT      = 0;
    xcmc = 0;
    
    %% Calculated Variables
    
    gamma = theta - alpha;
    qbar = 0.5*rho*V^2;
    
    % Thrust
    T = (550*th)*eta/max(V,1e-6);
    MT = eT*T;
    
    % Lift
    denL = 1 - qbar*S*(cbar/(2*V))*(CL_adot/(m*V));
    
    numL = qbar*S*( ...
        CL0 + CLa*alpha + CL_el*el ...
      + (cbar/(2*V))*CL_q*q ...
      + (cbar/(2*V))*CL_adot*( q + (W*cos(gamma) - T*sin(alpha - epsilon))/(m*V) ) );
    
    L = numL/denL;
    
    CL = L/(qbar*S);
    
    % Drag
    CD = CDM + k*(CL - CLDM)^2;
    D  = qbar*S*CD;
    
    % Moment
    alphadot = q + (-L + W*cos(gamma) - T*sin(alpha - epsilon))/(m*V);
    
    CM = CM_R0 + CM_Ra*alpha + CM_Rel*el ...
       + (cbar/(2*V))*CM_q*q ...
       + (cbar/(2*V))*CM_adot*alphadot ...
       - CL*xcmc;
    
    M = qbar*S*cbar*CM;
    
    %% State Derivatives
    Vdot     = (-D - W*sin(gamma) + T*cos(alpha - epsilon))/m;
    qdot     = (M + MT)/Jz;
    thetadot = q;
    %alphadot = q + (L + W*cos(gamma) - T*sin(alpha - epsilon))/(m*V);
    pdot     = V*cos(gamma);
    hdot     = V*sin(gamma);
    
    %% Outputs
    xdot = [Vdot; alphadot; qdot; thetadot; pdot; hdot];

end