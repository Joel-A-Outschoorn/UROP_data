warning ('off')
dJ = KS_MSS_prec(0)

function dJ = KS_MSS_prec(c)
T      = 10;
T_min  = -1000;
dt     = 1e-1;
t_span = (T_min:dt:T);
P      = 2;            % Number of segments
tp     = (0:T/P:T)';
L      = 128;          %   Length
N      = 127;          % # Interior nodes, 1,2,...,n
dx     = L/(N+1);
% init_cond = rand(N,1);
init_cond = linspace(1,1,N);
x      = dx:dx:L-dx;
%c      = 0;            % KS Parameter
maxit  = N*P;          % Max GMRES solver iterations

%%% Preconditioner parameters
No_sing = 15;               % Number of singular modes to compute
Tol     = 2;                
Maxit   = 1;                % Number of algorithm iterations
Subspace = No_sing+2;       % Subspace size
tol    = 1e-5;              % Schur Complement solver tolerance
gamma   = 0.05;              % Regularisation parameter

% t_span for the constraints and adjoints
for p=1:1:length(tp)-1
    tspan_con = tp(p):dt:tp(p+1);                                            
    t_span_con(:,p) = tspan_con';
end

opts   = odeset('RelTol',1e-4,'AbsTol',1e-4);
[t,u]  = ode45(@(t,u) ks_solve(t,u,dx,N,c), t_span, init_cond, opts);
t0_loc = ceil((length(t))*(-T_min/(T-T_min)));
t      = t(t0_loc:end);
u      = u(t0_loc:end,:);

% Compute f(u,t)
f1 = (2/(dx^2)-7/(dx^4))*u(:,1)  + (-(2*c+u(:,2))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,2) - (1/(dx^4))*u(:,3);
f2 = ((2*c+u(:,1))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,1) + (2/(dx^2)-6/(dx^4))*u(:,2) + (-(2*c+u(:,3))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,3) -(1/(dx^4))*u(:,4);
for i = 3:N-2
    f = (-1/(dx^4))*u(:,i-2) + ((2*c+u(:,i-1))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,i-1)  +  (2/(dx^2)-6/(dx^4))*u(:,i)  +  (-(2*c+u(:,i+1))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,i+1)  -  (1/(dx^4))*u(:,i+2);
    F(:,i-2) =f;
end
fNminus1 = -(1/(dx^4))*u(:,N-3) + ((2*c+u(:,N-2))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,N-2) +  (2/(dx^2)-6/(dx^4))*u(:,N-1)  +   (-(2*c+u(:,N))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,N);
fN = -(1/(dx^4))*u(:,N-2)   + ((2*c+u(:,N-1))/(4*dx)-1/(dx^2)+4/(dx^4)).*u(:,N-1)   +  (2/(dx^2)-7/(dx^4))*u(:,N);
F = [f1 f2 F fNminus1 fN];
% 

Init_cond = zeros(N,1);     % Zero initial condition to compute RHS
Prj = cell(P,1);
f_p = cell(P,1);

% Loop to compute RHS and projections
for p = 1:1:P
    F_p = F(ceil(1+(T/dt)*(p/P)),:)';     %f(u,t) at the checkpoints (p=1,2,...,P)
    
    [tc,v] = ode45(@(tc,v) con_solve(tc,v,t,u,dx,N,c), t_span_con(:,p), Init_cond, opts);   % Solve Constraint Equations with zero IC
    b = (eye(N) - (F_p*F_p')/(F_p'*F_p))*v(end,:)';                                         % Apply projection to zero IC terms
    C{p,:} = b;                                                                             % Store RHS
    Prj{p} = eye(N) - (F_p*F_p')/(F_p'*F_p);                                                % Store projections 
    f_p{p} = F_p;                                                                           % Store f(u,t) at the checkpoints (p=1,2,...,P)
    
    [U,S,V] = svds(@SVD_phi,[N N],No_sing,'largest','Tolerance',Tol,'MaxIterations',Maxit,'SubspaceDimension',Subspace); % Compute Singular modes
    M2 = U*(S^-2)*U' + (eye(N)-U*U');  % Preconditioner submatrices
    M_2{p} = M2;
end

Prj0 = mat2cell(eye(N) - (F(1,:)'*F(1,:))/(F(1,:)*F(1,:)'),N,N);            % Projection at t=0
Prj  = [Prj0;Prj];                                                          % Store projections at checkpoints (t0,t1,...,tP)
C    = -cell2mat(C);                                                        % RHS
M2   = blkdiag(M_2{1:end});                                                 % Form the full block diagonal preconditioner

J_x    = (1/L)*trapz(x,u,2);                                                % Space average             
J_xt   = (1/T)*trapz(t,J_x);                                                % Space-time average 

[ws,flag,relres1,iter,resvec1] = gmres(@Sw_PCvector,M2*C,[],tol,maxit);     % Solve Schur Complement using GMRES

% Compute  adjoint solution
for p =1:P
[ta,w] = ode45(@(ta,w) adj_solve(ta,w,t,u,dx,N,c), flipud(t_span_con(:,p)), Prj{p+1}*ws(1+N*(p-1):N*p), opts);   % Solve adjoint Equations 
w_end      = w(end,:)';
phiTw{p,:} = w_end;                                       % Adjoint Propagator MATVECs
end
phiT_w = cell2mat(phiTw);

vs    = -[-phiT_w(1:N); ws(1:N*(P-1))-phiT_w(N+1:end); ws(end-(N-1):end)];  % Constraint solution vector

% Compute derivative of time-average
for p = 0:1:P-1
    [tc,v_dash] = ode45(@(tc,v) con_solve(tc,v,t,u,dx,N,c),  t_span_con(:,p+1), vs(1+N*p:N*(p+1)), opts);
    
    f   = F(1-p+p*(length(t_span_con(:,1))):-p+(p+1)*(length(t_span_con(:,1))),:);
    vTf = dot(v_dash,f,2);
    f_norm2 = sum(f.*f,2);
    
    dJ_ds_p(p+1,:) = (1/T)*(sum(trapz(tc,(dx/L)*v_dash,1)) + ((f_p{p+1}'*v_dash(end,:)')/(f_p{p+1}'*f_p{p+1}))*(J_xt-J_x(ceil(1+(T/dt)*((p+1)/P)))));  % Sensitivity Integral
  
    v = v_dash -(vTf./f_norm2).*f;
end

dJ = sum(dJ_ds_p);   % Sensitivity (d<J_x>/dc)

% figure(3)
% semilogy(0:1:(length(resvec1)-1),resvec1,'-');
% hold on
% xlabel('Iteration number');
% ylabel('Residual');
% set(gca,'fontsize',13)

function MSw = eigPC(W)

for p =1:P
[ta,w] = ode45(@(ta,w) adj_solve(ta,w,t,u,dx,N,c), flipud(t_span_con(:,p)), Prj{p+1}*W(1+N*(p-1):N*p), opts);   % Solve adjoint Equations 
w_end      = w(end,:)';
phiTw{p,:} = w_end;                                       % Adjoint Propagator MATVECs
end
phiT_w = cell2mat(phiTw);
ATw    = [-phiT_w(1:N); W(1:N*(P-1))-phiT_w(N+1:end); W(end-(N-1):end)];  % MATVEC Product (A'w)

for p=1:P
[tc,v] = ode45(@(tc,v) con_hom_solve(tc,v,t,u,dx,N,c), t_span_con(:,p), ATw(1+N*(p-1):N*p), opts);   % Solve homogeneous Constraint Equations  
v_end   = v(end,:)';                                        % Constraint Propagator MATVECs
Sw{p,:} = -(Prj{p+1}*v_end)+ATw(1+(N*p):N*(p+1));           % MATVEC product (AA'w) 
end
Sw = cell2mat(Sw);
MSw = (gamma*eye(N*P))*W + M2*Sw;
end

function Ax = SVD_A(x,tflag)
         if strcmp(tflag,'notransp')
            for p=1:P
            [tc,v] = ode45(@(tc,v) con_hom_solve(tc,v,t,u,dx,N,c), t_span_con(:,p), x(1+N*(p-1):N*p), opts);   % Solve homogeneous Constraint Equations  
            v_end   = v(end,:)';                                        % Constraint Propagator MATVECs
            F_p = F(ceil(1+(T/dt)*(p/P)),:)';     %f(u,t) at the checkpoints (p=1,2,...,P)
            Ax{p,:} = -((eye(N) - (F_p*F_p')/(F_p'*F_p))*v_end)+x(1+(N*p):N*(p+1));           % MATVEC product (AA'w) 
            end
         Ax = cell2mat(Ax);                                      % Constraint Propagator MATVECs 
         else
         for p =1:P
              F_p = F(ceil(1+(T/dt)*(p/P)),:)';     %f(u,t) at the checkpoints (p=1,2,...,P) 
              [ta,w] = ode45(@(ta,w) adj_solve(ta,w,t,u,dx,N,c), flipud(t_span_con(:,p)), (eye(N) - (F_p*F_p')/(F_p'*F_p))*x(1+N*(p-1):N*p), opts);   % Solve adjoint Equations 
               w_end      = w(end,:)';
               phiTw{p,:} = w_end;                                       % Adjoint Propagator MATVECs
         end
               phiT_w = cell2mat(phiTw);
               Ax    = [-phiT_w(1:N); x(1:N*(P-1))-phiT_w(N+1:end); x(end-(N-1):end)];  % MATVEC Product (A'w)
         end
 end

function Phix = SVD_phi(x,tflag)
         if strcmp(tflag,'notransp')
         [tc,v] = ode45(@(tc,v) con_hom_solve(tc,v,t,u,dx,N,c), t_span_con(:,p), x, opts);   % Solve homogeneous Constraint Equations  
         Phix   = Prj{p}*v(end,:)';                                        % Constraint Propagator MATVECs 
         else
         [ta,w] = ode45(@(ta,w) adj_solve(ta,w,t,u,dx,N,c), flipud(t_span_con(:,p)), Prj{p}*x, opts);   % Solve adjoint Equations 
         Phix   = w(end,:)'; 
         end
 end

function Sw = Sw_PCvector(W)
    
    for p =1:P
    [ta,w] = ode45(@(ta,w) adj_solve(ta,w,t,u,dx,N,c), flipud(t_span_con(:,p)), Prj{p+1}*W(1+N*(p-1):N*p), opts);   % Solve adjoint Equations 
    w_end      = w(end,:)';
    phiTw{p,:} = w_end;                                       % Adjoint Propagator MATVECs
    end
    phiT_w = cell2mat(phiTw);
    ATw    = [-phiT_w(1:N); W(1:N*(P-1))-phiT_w(N+1:end); W(end-(N-1):end)];  % MATVEC Product (A'w)
    
    for p=1:P
    [tc,v] = ode45(@(tc,v) con_hom_solve(tc,v,t,u,dx,N,c), t_span_con(:,p), ATw(1+N*(p-1):N*p), opts);   % Solve homogeneous Constraint Equations  
    v_end   = v(end,:)';                                        % Constraint Propagator MATVECs
    Sw{p,:} = -(Prj{p+1}*v_end)+ATw(1+(N*p):N*(p+1));           % MATVEC product (AA'w) 
    end
    Sw = cell2mat(Sw);
    Sw = (gamma*eye(N*P))*W + M2*Sw;

end

function Sw = Sw_vector(W)
    
    for p =1:P
    [ta,w] = ode45(@(ta,w) adj_solve(ta,w,t,u,dx,N,c), flipud(t_span_con(:,p)), Prj{p+1}*W(1+N*(p-1):N*p), opts);   % Solve adjoint Equations 
    w_end      = w(end,:)';
    phiTw{p,:} = w_end;                                       % Adjoint Propagator MATVECs
    end
    phiT_w  = cell2mat(phiTw);
    ATw     = [-phiT_w(1:N); W(1:N*(P-1))-phiT_w(N+1:end); W(end-(N-1):end)];  % MATVEC Product (A'w)
    
    for p=1:P
    [tc,v]  = ode45(@(tc,v) con_hom_solve(tc,v,t,u,dx,N,c), t_span_con(:,p), ATw(1+N*(p-1):N*p), opts);   % Solve homogeneous Constraint Equations  
    v_end   = v(end,:)';                                        % Constraint Propagator MATVECs
    Sw{p,:} = -(Prj{p+1}*v_end)+ATw(1+(N*p):N*(p+1));           % MATVEC product (AA'w) 
    end
    Sw = cell2mat(Sw);
end

function dudt = ks_solve(t,u,dx,N,c)
       
         dudt = zeros(N,1);
         dudt(1) = (2/(dx^2)-7/(dx^4))*u(1)                   + (-(2*c+u(2))/(4*dx)-1/(dx^2)+4/(dx^4))*u(2) - (1/(dx^4))*u(3);
         dudt(2) = ((2*c+u(1))/(4*dx)-1/(dx^2)+4/(dx^4))*u(1) + (2/(dx^2)-6/(dx^4))*u(2) + (-(2*c+u(3))/(4*dx)-1/(dx^2)+4/(dx^4))*u(3) -(1/(dx^4))*u(4);
         for i = 3:N-2
         dudt(i) = (-1/(dx^4))*u(i-2) + ((2*c+u(i-1))/(4*dx)-1/(dx^2)+4/(dx^4))*u(i-1)  +  (2/(dx^2)-6/(dx^4))*u(i)  +  (-(2*c+u(i+1))/(4*dx)-1/(dx^2)+4/(dx^4))*u(i+1)  -  (1/(dx^4))*u(i+2); 
         end
         dudt(N-1) = -(1/(dx^4))*u(N-3) + ((2*c+u(N-2))/(4*dx)-1/(dx^2)+4/(dx^4))*u(N-2) +  (2/(dx^2)-6/(dx^4))*u(N-1)  +   (-(2*c+u(N))/(4*dx)-1/(dx^2)+4/(dx^4))*u(N);
         dudt(N) = -(1/(dx^4))*u(N-2)   + ((2*c+u(N-1))/(4*dx)-1/(dx^2)+4/(dx^4))*u(N-1)   +  (2/(dx^2)-7/(dx^4))*u(N);
end

function dvdt = con_solve(tc,v,t,u,dx,N,c)
         
         u    = interp1(t,u,tc);
         
         dvdt    = zeros(N,1);
         dvdt(1) = (2/(dx^2)-7/(dx^4))*v(1) + (-(c+u(2))/(2*dx)-1/(dx^2)+4/(dx^4))*v(2) - (1/(dx^4))*v(3) - u(2)/(2*dx);
         dvdt(2) = ((c+u(1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(1) + (2/(dx^2)-6/(dx^4))*v(2) + (-(c+u(3))/(2*dx)-1/(dx^2)+4/(dx^4))*v(3) -(1/(dx^4))*v(4) + u(1)/(2*dx) - u(3)/(2*dx);
         for i = 3:N-2
         dvdt(i) = (-1/(dx^4))*v(i-2) + ((c+u(i-1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(i-1)  +  (2/(dx^2)-6/(dx^4))*v(i)  +  (-(c+u(i+1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(i+1)  -  (1/(dx^4))*v(i+2) + u(i-1)/(2*dx) - u(i+1)/(2*dx); 
         end
         dvdt(N-1) = -(1/(dx^4))*v(N-3) + ((c+u(N-2))/(2*dx)-1/(dx^2)+4/(dx^4))*v(N-2) +  (2/(dx^2)-6/(dx^4))*v(N-1)  +   (-(c+u(N))/(2*dx)-1/(dx^2)+4/(dx^4))*v(N) + u(N-2)/(2*dx) - u(N)/(2*dx);
         dvdt(N) = -(1/(dx^4))*v(N-2)   + ((c+u(N-1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(N-1)   +  (2/(dx^2)-7/(dx^4))*v(N) + u(N-1)/(2*dx);
end

function dvdt = con_hom_solve(tc,v,t,u,dx,N,c)
         
         u    = interp1(t,u,tc);
         
         dvdt    = zeros(N,1);
         dvdt(1) = (2/(dx^2)-7/(dx^4))*v(1) + (-(c+u(2))/(2*dx)-1/(dx^2)+4/(dx^4))*v(2) - (1/(dx^4))*v(3);
         dvdt(2) = ((c+u(1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(1) + (2/(dx^2)-6/(dx^4))*v(2) + (-(c+u(3))/(2*dx)-1/(dx^2)+4/(dx^4))*v(3) -(1/(dx^4))*v(4);
         for i = 3:N-2
         dvdt(i) = (-1/(dx^4))*v(i-2) + ((c+u(i-1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(i-1)  +  (2/(dx^2)-6/(dx^4))*v(i)  +  (-(c+u(i+1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(i+1)  -  (1/(dx^4))*v(i+2); 
         end
         dvdt(N-1) = -(1/(dx^4))*v(N-3) + ((c+u(N-2))/(2*dx)-1/(dx^2)+4/(dx^4))*v(N-2) +  (2/(dx^2)-6/(dx^4))*v(N-1)  +   (-(c+u(N))/(2*dx)-1/(dx^2)+4/(dx^4))*v(N);
         dvdt(N) = -(1/(dx^4))*v(N-2)   + ((c+u(N-1))/(2*dx)-1/(dx^2)+4/(dx^4))*v(N-1)   +  (2/(dx^2)-7/(dx^4))*v(N);
end

function dwdt = adj_solve(ta,w,t,u,dx,N,c)
         
         u    = interp1(t,u,ta);
    
         dwdt    = zeros(N,1);
         dwdt(1) = -(2/(dx^2)-7/(dx^4))*w(1) - ((c+u(1))/(2*dx) - 1/(dx^2)+4/(dx^4))*w(2) + (1/(dx^4))*w(3);
         dwdt(2) = -(-(c+u(2))/(2*dx)-1/(dx^2)+4/(dx^4))*w(1) - (2/(dx^2)-6/(dx^4))*w(2) - ((c+u(2))/(2*dx)-1/(dx^2)+4/(dx^4))*w(3) + (1/(dx^4))*w(4);
%          dwdt(3) = (1/(dx^4))*w(1) - (-(c+u(3))/(2*dx)-1/(dx^2)+4/(dx^4))*w(2) - (2/(dx^2)-6/(dx^4))*w(3) - ((c+u(3))/(2*dx)-1/(dx^2)+4/(dx^4))*w(4) + (1/(dx^4))*w(5);
%          dwdt(4) = (1/(dx^4))*w(2) - (-(c+u(4))/(2*dx)-1/(dx^2)+4/(dx^4))*w(3) - (2/(dx^2)-6/(dx^4))*w(4) - ((c+u(4))/(2*dx)-1/(dx^2)+4/(dx^4))*w(5) + (1/(dx^4))*w(6);
         for i   = 3:N-2
         dwdt(i) =  (1/(dx^4))*w(i-2) - (-(c+u(i))/(2*dx)-1/(dx^2)+4/(dx^4))*w(i-1) - (2/(dx^2)-6/(dx^4))*w(i) - ((c+u(i))/(2*dx)-1/(dx^2)+4/(dx^4))*w(i+1) + (1/(dx^4))*w(i+2);
         end
%          dwdt(N-3) = (1/(dx^4))*w(N-5) - (-(c+u(N-3))/(2*dx)-1/(dx^2)+4/(dx^4))*w(N-4) - (2/(dx^2)-6/(dx^4))*w(N-3) - ((c+u(N-3))/(2*dx)-1/(dx^2)+4/(dx^4))*w(N-2) + (1/(dx^4))*w(N-1);
%          dwdt(N-2) = (1/(dx^4))*w(N-4) - (-(c+u(N-2))/(2*dx)-1/(dx^2)+4/(dx^4))*w(N-3) - (2/(dx^2)-6/(dx^4))*w(N-2) - ((c+u(N-2))/(2*dx)-1/(dx^2)+4/(dx^4))*w(N-1) + (1/(dx^4))*w(N);
         dwdt(N-1) = (1/(dx^4))*w(N-3) - (-(c+u(N-1))/(2*dx)-1/(dx^2)+4/(dx^4))*w(N-2) - (2/(dx^2)-6/(dx^4))*w(N-1) - ((c+u(N-1))/(2*dx)-1/(dx^2)+4/(dx^4))*w(N);
         dwdt(N)   = (1/(dx^4))*w(N-2) - (-(c+u(N))/(2*dx)-1/(dx^2)+4/(dx^4))*w(N-1) - (2/(dx^2)-7/(dx^4))*w(N);
end

end