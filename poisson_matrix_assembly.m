function [A,b] = poisson_matrix_assembly(nx,dx,f,BC_type,BC_value)
%poisson_matrix_assembly
%   this function assembles the linear equation system matrix
%   for the 1D Poisson equation
A = -2*eye(nx)/dx/dx;
A(2:end,2:end) = A(2:end,2:end) + diag(ones(nx-2,1)/dx/dx,1);
A(1:end-1,1:end-1) = A(1:end-1,1:end-1) + diag(ones(nx-2,1)/dx/dx,-1);

% source vector
b = ones(nx,1);
b(2:end-1) = f;

% left boundary condition treatment
if BC_type{1} == 'DBC'
    A(1,1) = 1;
    b(1) = BC_value(1);
elseif BC_type{1} == 'NBC'
    A(1,1) = 1/dx;
    A(1,2) = -1/dx;
    b(1) = BC_value(1);
else
    error(['unrecognised boundary type: ',BC_type{1}])
end

% right boundary condition treatment
if BC_type{2} == 'DBC'
    A(end,end) = 1;
    b(end) = BC_value(2);
elseif BC_type{2} == 'NBC'
    A(end,end) = -1/dx;
    A(end,end-1) = 1/dx;
    b(end) = BC_value(2);
else
    error(['unrecognised boundary type: ',BC_type{2}])
end
end