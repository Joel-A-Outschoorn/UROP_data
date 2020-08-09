clear
clc
close all

sigma=10; rho=28; eps = 0.8;
beta_min = 2; beta_max = 2.8; N = 10;
beta = 8/3;
%ii=1; dJ = zeros(N,2);; 
a1 = normrnd(1,0.3); a2 = normrnd(1,0.01); a3 = normrnd(1,0.01);
init_cond = [a1 a2 a3];

% 3D plot
T = 1000;
[J_mean,x,y,z] = Lorenzz(sigma,rho,beta,init_cond,T);
figure
plot3(x,y,z);

%Plot of convergence of J
T_plot = 200:5000:50200;
J_mean_plot1 = zeros;
for i = 1:length(T_plot)
    J_mean_plot1(i) = Lorenzz(sigma,rho,beta,init_cond,T_plot(i));
end
figure
plot(T_plot,J_mean_plot1)

%Plot of change in J with beta
beta_plot = linspace(beta_min,beta_max,N);
T_plot2 = [100 500 1000];
J_mean_plot2 = zeros;
for i = 1:length(T_plot2)
    for j = 1:length(beta_plot)
        J_mean_plot2(i,j) = Lorenzz(sigma,rho,beta_plot(j),init_cond,T_plot2(i));
    end
end
figure
plot(beta_plot,J_mean_plot2(1,:),'r')
hold on
plot(beta_plot,J_mean_plot2(2,:),'b')
plot(beta_plot,J_mean_plot2(3,:),'g')

return

function dJ = findif(sigma,rho,beta,eps)
     Jip1  = Lorenzz(sigma,rho,beta+eps);
     Jim1  = Lorenzz(sigma,rho,beta-eps);
     dJ=Jip1-Jim1; dJ=dJ/eps; dJ=dJ/2;
     %Jip2 = Lorenzz(sigma,rho+2*eps,beta);
     %Jip1 = Lorenzz(sigma,rho+  eps,beta);
     %Jim1 = Lorenzz(sigma,rho-  eps,beta);
     %Jim2 = Lorenzz(sigma,rho-2*eps,beta);
     %dJ =-Jip2+8*Jip1-8*Jim1+Jim2; dJ = dJ/eps; dJ = dJ/12;
end


function [J_mean,x,y,z] = Lorenzz(sigma,rho,beta,init_cond,T)
dt    = 1e-2;                                 % Integrator interpolation time step
T_min = -50;                                  % Integrate Lorenz from Tmin to T. Discard transient section (Tmin to 0).
P         = 1000;                                % Number of segments
tp        = (0:T/P:T)';
tspan     = (T_min:dt:T);                     % Time interval for the Lorenz
%init_cond = ones(3,1);                        % Lorenz Initial Condition Vector
%beta  = 8/3;
opts  = odeset('RelTol',1e-5,'AbsTol',1e-5);  % Integrator tolerances
for p=1:1:length(tp)-1
    tspan_con = tp(p):dt:tp(p+1);                                            
    t_span_con(:,p) = tspan_con';
end
[t,u] = ode45(@(t,u) lorenz_solve(t,u,sigma,rho,beta), tspan, init_cond, opts);  % Solve the Lorenz System
% Discard transient section and compute f1, f2, f3
t0_loc  = ceil((length(t))*(-T_min/(T-T_min)));
t       = t(t0_loc:end);
x       = u(t0_loc:end,1); y  = u(t0_loc:end,2); z  = u(t0_loc:end,3);           % Lorenz Solutions
f1      = sigma*(y-x);     f2 = x.*(rho-z)-y;    f3 = x.*y - beta*z;             % f(u,t)
df2_ds    = x;                                                                   % Forcing term df/ds
J         = z;                                                                   % Instantaneous cost function
J_mean    = (1/T)*simpsons(J,0,T,[]);   
% Time averaged cost function
end

% Lorenz solver. typical values: rho = 28; sigma = 10; beta = 8/3;
function dudt = lorenz_solve(t,u,sigma,rho,beta)
dudt = zeros(3,1);
dudt(1) = sigma*(u(2) - u(1));
dudt(2) = u(1)*(rho - u(3)) - u(2);
dudt(3) = u(1)*u(2) - beta*u(3);
end