function [x, flag, relres, iter, resvec] = jacobi_solver(A, b, tol, maxit, x0)
%JACOBI_SOLVER Solve Ax=b using the Jacobi iterative method.
%
% Inputs:
%   A      - (n x n) matrix
%   b      - (n x 1) right-hand side
%   tol    - stopping tolerance on relative residual (default 1e-8)
%   maxit  - maximum iterations (default 500)
%   x0     - initial guess (default zeros(n,1))
%
% Outputs:
%   x      - approximate solution
%   flag   - 0 if converged, 1 if not converged, 2 if breakdown/invalid
%   relres - final relative residual norm(b-Ax)/norm(b)
%   iter   - iterations performed
%   resvec - residual history (relative)

    if nargin < 3 || isempty(tol),   tol = 1e-8; end
    if nargin < 4 || isempty(maxit), maxit = 500; end
    n = length(b);
    if nargin < 5 || isempty(x0),    x0 = zeros(n,1); end

    % Basic checks
    flag = 1;
    iter = 0;
    x = x0;

    if ~ismatrix(A) || size(A,1) ~= size(A,2) || size(A,1) ~= n
        flag = 2;
        error('A must be square and compatible with b.');
    end

    D = diag(A);
    if any(D == 0)
        flag = 2;
        error('Jacobi breakdown: A has a zero diagonal entry.');
    end

    % Decompose: A = D + R, where D is diagonal and R is the remainder
    R = A - diag(D);

    bnorm = norm(b);
    if bnorm == 0
        bnorm = 1; % avoid divide-by-zero in relative residual
    end

    resvec = zeros(maxit+1,1);
    r = b - A*x;
    resvec(1) = norm(r) / bnorm;

    if resvec(1) <= tol
        flag = 0;
        relres = resvec(1);
        resvec = resvec(1);
        return;
    end

    for k = 1:maxit
        % Jacobi update: x^{k+1} = D^{-1}(b - R x^k)
        x_new = (b - R*x) ./ D;

        % Residual and stopping test
        r = b - A*x_new;
        resvec(k+1) = norm(r) / bnorm;

        x = x_new;
        iter = k;

        if resvec(k+1) <= tol
            flag = 0;
            break;
        end
    end

    relres = resvec(iter+1);
    resvec = resvec(1:iter+1);
end
