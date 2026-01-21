function [u] = linear_system_solver4Poisson(A,b,solver_type,niter)
%linear_system_solver4Poisson
%   provides linear system solver
if strcmp(solver_type, 'direct')
    u = A\b;
elseif strcmp(solver_type, 'indirect')
    [u, flag, relres, iter, resvec] = jacobi_solver(A, b, 1e-10, niter);
else
    error(['unrecognised solver type: ',solver_type])
end

end