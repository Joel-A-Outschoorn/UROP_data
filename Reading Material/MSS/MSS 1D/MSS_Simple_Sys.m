% SENSITIVITY ANALYSIS OF SIMPLE ODE SYSTEM
% dydt = -k*y
% J = y
% J_mean = (1/T) * int_0^T(y)dt
% Design variable,s -> s = k

clear
clc
close all

% Constants
k = 0.25;
T = 30;
P = T;
IC = 0; %Initial condition
opts = odeset('RelTol',1e-5,'AbsTol',1e-5); %Integrator tolerances

%% Analytic solution 
t = 0:0.001:T;
y = exp(-k*t);
J = y;
J_bar = (1/T)*simpsons(J,0,T,30);
f = -k*y;

plot(t,y,'-r')

dJ_analytic = (1/k)*exp(-k*T) - (1/(k^2*T))*(1 - exp(-k*T));

%% MSS 
dJ_MSS = MSS(10,0.8);

function dJ_MSS = MSS(k,e)
% Constants
T = 30; P = T; IC = 0; %Initial condition
dt = 1e-3; %Time step
tol = 1e-10; %GMRES tol
maxit = P; %No. of GMRES iterations
tp = (0:T/P:T)';
opts = odeset('RelTol',1e-5,'AbsTol',1e-5); %Integrator tolerances

    %Time span for constraints and adjoint
    for p = 1:length(tp)-1
        tspan_constraint = tp(p):dt:tp(p+1);
        t_span_constraint(:,p) = tspan_constraint';
    end
    
%Solving ODE system
t = 0:0.001:T;
y = exp(-k*t);
%f(u,t)
f = -k*y;

J = y;
J_bar = (1/T)*simpsons(J,0,T,30);
IC_MSS = 0;

%Initialisation
Prj = cell(P,1);
PHI = cell(P,1);
F = cell(P,1);
f_p = cell(P,1);
c = cell(P,1);

    %Loop to compute projections at each point
    for p =1:P
        %Finding f at each point,p
        F1 = f(find(t == p));
        
        %Solving constraint equation with zero IC
        [tc,v] = ode45(@(tc,v) constraint_solve(tc,v,t,y,k,IC_MSS),t_span_constraint(:,p),IC_MSS,opts);
        %v at point,p
        b = v(end,:)';
        %Solving constraint equation with slight perturbation of IC
        [tc,v] = ode45(@(tc,v) constraint_solve(tc,v,t,y,k,IC_MSS),t_span_constraint(:,p),e,opts);
        %v at point,p
        phi1 = (v(end,:)'-b)/e;     
        
        v_end = (eye(1) - (F1*F1')/(F1'*F1))*b;
        phi = (eye(1) - (F1*F1')/(F1'*F1))*phi1;
        
        %Storing projections
        PHI{p} = phi;
        c{p,:} = v_end;
        Prj{p} = (eye(1) - (F1*F1')/(F1'*F1)); 
        F{p} = F1;
        f_p{p} = F;
    end

%Constraint matrix 
A = [-blkdiag(PHI{1:end}) zeros(1*P,1)]; A = A + [zeros(1*P,1) eye(1*P)];  

%Schur component
Sch = A*A';

%Projection at t=0
Prj0 = mat2cell(eye(1) - (f(1)*f(1))/(f(1)*f(1)),1,1);     % Projection at t=0
Prj  = [Prj0;Prj];                                        % Store projections at checkpoints t0, t1,...,tP

c = cell2mat(c);

%Solving Schur complement using GMRES
[L,flag,relres,iter,resvec] = gmres(Sch,c,[],tol,maxit);

%Constraint solution vector
Xs = A'*L; 

    %Final solution
    for p = 0:P-1
        %Solving constraint equation with solved IC
        [tc,v] = ode45(@(tc,v) constraint_solve(tc,v,t,y,k,Xs(p+1)),t_span_constraint(:,p+1),Xs(p+1),opts);
        v_dash = v;
        
        %Sensitivity integral
        % dJ_p(p+1,:) = (1/T)*simpsons(v_dash,(p/P)*T,((p+1)/P)*T,[]) + (1/T)*((f_p{p+1}'*v(end,:)')/(f_p{p+1}'*f_p{p+1}))*(J_bar - J(find(t == (p+1))))
        dJ_p(p+1,:) = (1/T)*simpsons(v_dash,(p/P)*T,((p+1)/P)*T,[]) + (1/T)*((sum(cell2mat(f_p{p+1})'*v(end,:)))/sum(cell2mat(f_p{p+1})'*cell2mat(f_p{p+1})))*(J_bar - J(find(t == (p+1))));
    end
    
%FINAL VALUE
dJ_MSS = sum(dJ_p);

    % Constraint equation
    function dvdt = constraint_solve(tc,v,t,y,k,IC)

        y = interp1(t,y,tc);

        dvdt = IC;
        dvdt = -k*v - y;
    end 
end