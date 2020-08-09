
tic

% format long;
% clear
% clc
close all

% h=0.0001; T = 0.5; t=0:h:T; %time step and vector containint
% id = 100; ytarget = 10; %id =  number of design variables
% kvector = 3*ones(id,1);
% y0vector = ones(id,1); %here we have 10 kappas. each can be different but here they are all equal to each other
[J,y] = primal(h,t,T,kvector,y0vector,ytarget); %here we solve the primal and return the field y (which is needed for the adjoint) and the objective J
dJ = adjoint(h,t,T,kvector,y,ytarget); %here dJ is a vector containing the sensitivities of J w.r.t. the vector kappa

figure
for i = 1:id
plot(t,y(:,i));
hold on
end

toc

function dF = adjoint(h,t,T,kvector,y,ytarget)
%
f=@(ay,term,kappa) kappa*ay + term;  
Nt = length(t); 
N = length(kvector); 
ay = zeros(Nt,N); % Impose Adjoint Boundary Conditions & preallocation
for j = 1:N
    kap = kvector(j);
    for ii=1:(Nt-1) %loop
        i = Nt - ii;   
       sumsum = (y(i,j)-ytarget)/N;
%      RK4 coefficients
        k1=f(ay(ii),sumsum,kap);
%         k2=f(ay(ii)+(0.5*k1*h),sumsum,kap) 
%         k3=f(ay(ii)+(0.5*k2*h),sumsum,kap)
%         k4=f(ay(ii)+k3*h,sumsum,kap)
%         pause
%     Integrate for ay
        ay(ii+1,j) = ay(ii,j) - h*(k1);% + 2*k2 + 2*k3 + k4)/6;
    end
end
dF = zeros(N,1);
for iii = 1:N
    tempp = zeros(Nt,1);
    for i = 1:Nt
        kk = Nt - i + 1;
        tempp(i) = ay(i,iii)*y(kk,iii);
    end
    dF(iii) = simpsons(tempp,0,T,[]);
end
plot(t,ay)
end 


function [J,y] = primal(h,t,T,kvector,y0vector,ytarget)
%
%kvevtor contains N kappa values, y0 vector contains N initial conditions
%
N = length(kvector); Nt = length(t);
y = zeros(Nt,N);
%
%solve the primal at each dimension
for ii = 1:N     
    y(:,ii) = rkk4(kvector(ii),h,t,y0vector(ii));
end
%
%find objective function
%
sum = 0;
for i = 1:N
    temp = zeros(Nt,1);
    for j = 1:Nt
        temp(j) = 0.5*(y(j,i)-ytarget)^2;  
    end
    int = simpsons(temp,0,T,[]);
    sum = sum + int;   
end
J = sum/N;
end

function y = rkk4(kappa,h,t,y0) %with RK4
f=@(y) -kappa*y;  
y = zeros(length(t),1); 
y(1) = y0;
for i=1:(length(t)-1) %loop
    k1=f(y(i));
    k2=f(y(i)+(0.5*k1*h));
    k3=f(y(i)+(0.5*k2*h));
    k4=f(y(i)+k3*h);
    y(i+1) = y(i) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
end
end