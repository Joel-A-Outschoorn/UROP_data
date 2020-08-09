%Finite difference method
tic

clear 
clc
close all


h=0.0001; T = 0.5; t=0:h:T; % time step and vector containint
id = 100; ytarget = 10; % id =  number of design variables
kvector = 3*ones(id,1);  
% kvector = randi([1 10],id,1);  
y0vector = ones(id,1); % here we have 10 kappas. each can be different but here they are all equal to each other
% y0vector = randi([1 10],id,1); 
k = kvector;
J_k = primal(h,t,T,k,y0vector,ytarget); % here we solve the primal and return the field y (which is needed for the adjoint) and the objective J

%Perturbing k vector
dk = 0.005;
dJ_FD = zeros;
J_k_Plus_dk = zeros;
for i = 1:length(kvector)
    kvector_perturb = kvector;
    kvector_perturb(i) = dk + kvector(i);
    k_Plus_dk = kvector_perturb;
    J_k_Plus_dk(i) = primal(h,t,T,k_Plus_dk,y0vector,ytarget);
    
    %FD approximation of sensitivities
    dJ_FD(i,1) = (J_k_Plus_dk(i) - J_k)/dk;
end

toc

%Comparison
run('adjointFD.m')
Percent_diff = abs((dJ_FD - dJ)./ dJ) * 100;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONS 

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
