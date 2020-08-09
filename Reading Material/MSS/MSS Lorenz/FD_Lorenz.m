close all
clear
clc

tic
J_beta = MSS_Lorenzz(10,28,2.6,0.8)
db = 100;
beta_plus_db = 2.6 + db;
J_beta_plus_db = MSS_Lorenzz(10,28,beta_plus_db,0.8)

dJ = (J_beta_plus_db - beta_plus_db) / db;
toc

function J_beta = MSS_Lorenzz(sigma,rho,beta,e)

T = 30;  P = T;                         % Time Horizon
dt  = 1e-3;                                 % Time Step
T_min = -50;                                  % Integrate Lorenz from Tmin to T. Discard transient section.                               % Number of segments
tp        = (0:T/P:T)';
tspan     = (T_min:dt:T);                     % Time interval for the Lorenz
init_cond = ones(3,1);                        % Lorenz Initial Condition Vector
opts  = odeset('RelTol',1e-5,'AbsTol',1e-5);  % Integrator tolerances
tol   = 1e-10;                                % GMRES tol
maxit = 3*P;                                  % Number of GMRES Iterations
%e     = 1e-4;                                % Finite differ

% t_span for the constraints and adjoints
for p=1:1:length(tp)-1
    tspan_con = tp(p):dt:tp(p+1);                                            
    t_span_con(:,p) = tspan_con';
end

[t,u] = ode45(@(t,u) lorenz_solve(t,u,sigma,rho,beta), tspan, init_cond, opts);  % Solve the Lorenz

t0_loc  = ceil((length(t))*(-T_min/(T-T_min)));
t       = t(t0_loc:end);
x       = u(t0_loc:end,1); y  = u(t0_loc:end,2); z  = u(t0_loc:end,3);           % Lorenz Solutions
f1      = sigma*(y-x);     f2 = x.*(rho-z)-y;    f3 = x.*y - beta*z;             % f(u,t)

df3_ds    = -z;                                                                   % Forcing term df/ds
J         = z;
J_mean    = (1/T)*simpsons(J,0,T,[]);  
J_beta = J_mean;

function dudt = lorenz_solve(t,u,sigma,rho,beta)
dudt = zeros(3,1);
dudt(1) = sigma*(u(2) - u(1));
dudt(2) = u(1)*(rho - u(3)) - u(2);
dudt(3) = u(1)*u(2) - beta*u(3);
end

end
