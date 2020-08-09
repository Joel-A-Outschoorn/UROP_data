close all
clear
clc

tic
dJ = MSS_Lorenzz(10,28,2.6,0.8)
toc

function dJ = MSS_Lorenzz(sigma,rho,beta,e)
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
Init_cond = zeros(3,1);

Prj = cell(P,1);
PHI = cell(P,1);
F_1 = cell(P,1); F_2 = cell(P,1); F_3 = cell(P,1);
f = cell(P,1);
c = cell(P,1);

% Loop to compute RHS and projections
for p = 1:1:P
    % Find f at checkpoints p=1 to p=end
    F1 = f1(ceil(1+(T/dt)*(p/P)));  F2 = f2(ceil(1+(T/dt)*(p/P)));  F3 = f3(ceil(1+(T/dt)*(p/P)));
    F  = [F1 F2 F3]';
    
    %%% 
    [tc,v] = ode45(@(tc,v) con_solve(tc,v,t,x,y,z,rho,sigma,beta,df3_ds,Init_cond), t_span_con(:,p), Init_cond, opts);   % Solve Constraint Equations with zero IC
    b = v(end,:)';
    [tc,v] = ode45(@(tc,v) con_solve(tc,v,t,x,y,z,rho,sigma,beta,df3_ds,[e;0;0]), t_span_con(:,p), [e;0;0], opts);   % Solve Constraint Equations with zero IC
    phi1 = (v(end,:)'- b)/e;
    [tc,v] = ode45(@(tc,v) con_solve(tc,v,t,x,y,z,rho,sigma,beta,df3_ds,[0;e;0]), t_span_con(:,p), [0;e;0], opts);   % Solve Constraint Equations with zero IC
    phi2 = (v(end,:)'- b)/e;
    [tc,v] = ode45(@(tc,v) con_solve(tc,v,t,x,y,z,rho,sigma,beta,df3_ds,[0;0;e]), t_span_con(:,p), [0;0;e], opts);   % Solve Constraint Equations with zero IC
    phi3 = (v(end,:)'- b)/e;
    v_end  = (eye(3) - (F*F')/(F'*F))*b;                                                                     % Zero IC terms with projection applied
    phi = (eye(3) - (F*F')/(F'*F))*[phi1 phi2 phi3];

    PHI{p} = phi;
    c{p,:} = v_end;  

    Prj{p} = eye(3) - (F*F')/(F'*F);              % Store projections 
    F_1{p} = F1; F_2{p} = F2; F_3{p} = F3;
    f{p}   = F;
end

A = [-blkdiag(PHI{1:end}) zeros(3*P,3)]; A = A + [zeros(3*P,3) eye(3*P)];  % MSS Constraint Matrix

Sch = A*A';     % Schur Complement

Prj0 = mat2cell(eye(3) - ([f1(1);f2(1);f3(1)]*[f1(1) f2(1) f3(1)])/([f1(1) f2(1) f3(1)]*[f1(1);f2(1);f3(1)]),3,3);     % Projection at t=0
Prj  = [Prj0;Prj];                                        % Store projections at checkpoints t0, t1,...,tP

c   = cell2mat(c);

[L,flag,relres,iter,resvec] = gmres(Sch,c,[],tol,maxit);  % Solve schur complement using GMRES


Xs = A'*L;                                                % Constraint Solution Vector

% figure
% semilogy(0:1:(length(resvec)-1),resvec,'-'); xlabel('Iteration number');
% ylabel('Relative Residual');
% set(gca,'fontsize',13)
% hold on

% Compute Solution and find derivatives
for p = 0:1:P-1
    [tc,v] = ode45(@(tc,v) con_solve(tc,v,t,x,y,z,rho,sigma,beta,df3_ds,Xs(1+3*(p):3*(p+1))), t_span_con(:,p+1), Xs(1+3*(p):3*(p+1)), opts);   % Constraint Equations
    v1dash = v(:,1); v2dash = v(:,2); v3dash = v(:,3);
    
    dJ_ds_p(p+1,:) = (1/T)*(simpsons(v3dash,(p/P)*T,((p+1)/P)*T,[]) + ((f{p+1}'*v(end,:)')/(f{p+1}'*f{p+1}))*(J_mean-J(ceil(1+(T/dt)*((p+1)/P)))));  % Sensitivity Integral

    f_1 = f1(1-p+p*(length(t_span_con(:,1))):-p+(p+1)*(length(t_span_con(:,1))));
    f_2 = f2(1-p+p*(length(t_span_con(:,1))):-p+(p+1)*(length(t_span_con(:,1))));
    f_3 = f3(1-p+p*(length(t_span_con(:,1))):-p+(p+1)*(length(t_span_con(:,1))));

    fTvdash = v1dash.*f_1 + v2dash.*f_2 + v3dash.*f_3;
    f_norm2 = f_1.^2 + f_2.^2 + f_3.^2;
    
    v1 = v1dash - (fTvdash./f_norm2).*f_1;
    v2 = v2dash - (fTvdash./f_norm2).*f_2;
    v3 = v3dash - (fTvdash./f_norm2).*f_3;
    
%     figure
%     plot(tc,v1,'r')
%     hold on
%     plot(tc,v2,'b')
%     plot(tc,v3,'k')
%     title('Constraints')
end

dJ_ds = sum(dJ_ds_p);% Sensitivity (d<J>/dp)
dJ = dJ_ds;
% typical values: rho = 28; sigma = 10; beta = 8/3;
function dudt = lorenz_solve(t,u,sigma,rho,beta)
dudt = zeros(3,1);
dudt(1) = sigma*(u(2) - u(1));
dudt(2) = u(1)*(rho - u(3)) - u(2);
dudt(3) = u(1)*u(2) - beta*u(3);
end

% Function to solve Constraint Equations
function  dvdt = con_solve(tc,v,t,x,y,z,rho,sigma,beta,df3_ds,Init_cond)  
    
          x  = interp1(t,x,tc);
          y  = interp1(t,y,tc);
          z  = interp1(t,z,tc);
          df3_ds = interp1(t,df3_ds,tc);

          dvdt    = Init_cond;
          dvdt(1) = -sigma*v(1) + sigma*v(2);
          dvdt(2) = (rho-z)*v(1) -       v(2) -     x*v(3);
          dvdt(3) =       y*v(1) +     x*v(2) -  beta*v(3)  + df3_ds;
end

end