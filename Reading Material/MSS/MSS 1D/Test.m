% Constants
T = 30; P = T; IC = 0; %Initial condition
dt = 1e-3; %Time step
tol = 1e-10; %GMRES tol
maxit = 3*P; %No. of GMRES iterations
tp = (0:T/P:T)';

    %Time span for constraints and adjoint
    for p = 1:length(tp)-1
        tspan_constraint = tp(p):dt:tp(p+1);
        t_span_constraint(:,p) = tspan_constraint';
    end
    
A = [2 3; 1 2];
B = [4 4; 4 4]; x{1} = A; x{2} = B;

[blkdiag(x{1:end})]