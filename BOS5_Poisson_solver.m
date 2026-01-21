clear all;
close all;
clc;

%% physical parameters
% domain specification
x0 = 0;
x1 = 1;

% boundary conditions
% Dirichlet - DBC; Neumann - NBC;
BC_type = {'DBC','DBC'};
% unknown u for DBC, unknown derivative for NBC
BC_value = [0,0];


%% numerical configuration
% number of nodes
nx = 4;

solver_type = 'indirect';
niter = 200000;

x = linspace(x0,x1,nx);
dx = abs(x(2)-x(1));

% source term
f = 5;%sin(x(2:end-1));

%% matrix assembly
[A,b] = poisson_matrix_assembly(nx,dx,f,BC_type,BC_value);

%% obtaining numerical solution
tic;
u = linear_system_solver4Poisson(A,b,solver_type,niter);
toc

%% obtaining analytical solution
a = f/2;
C2 = BC_value(1);
C1 = BC_value(2) - f/2 - C2;

x_analyt = x.';%linspace(x0,x1,1001);
u_analyt = a*x_analyt.^2 + C1*x_analyt + C2;

err = abs(u - u_analyt);
disp(max(err));

%% visualisation
% --- Figure: 1 row x 2 columns ---
fig = figure('Units','centimeters','Position',[5 5 17 5]);   % center-ish, 17×5 cm

set(fig,'PaperUnits','centimeters')
set(fig,'PaperSize',[17 5])
set(fig,'PaperPosition',[0 0 17 5])
set(fig,'PaperPositionMode','manual')


% Global style controls
fs = 10; % font size (change as you like)

tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- Left plot: solution vs analytic ---
nexttile
plot(x, u, 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k'); hold on
plot(x_analyt, u_analyt, '-r', 'LineWidth', 2);
grid on
xlabel('$x$', 'Interpreter','latex', 'FontSize', fs)
ylabel('$u$', 'Interpreter','latex', 'FontSize', fs)
legend({'Numerical','Analytical'}, 'Interpreter','latex', 'FontSize', fs, 'Location','north')
set(gca,'FontSize',fs,'TickLabelInterpreter','latex')

annotation('textbox',[0.1 0.4 0.05 0.05], ...
    'String','(a)','Interpreter','latex', ...
    'EdgeColor','none','FontSize',12,'FontWeight','bold');

% --- Right plot: pointwise L2 error (actually absolute error) ---
nexttile
plot(x, err, 'LineWidth', 2);
grid on
xlabel('$x$', 'Interpreter','latex', 'FontSize', fs)
ylabel('$|u-u_{\mathrm{analyt}}|$', 'Interpreter','latex', 'FontSize', fs)
legend({'Error'}, 'Interpreter','latex', 'FontSize', fs, 'Location','south')
set(gca,'FontSize',fs,'TickLabelInterpreter','latex')

annotation('textbox',[0.6 0.4 0.05 0.05], ...
    'String','(b)','Interpreter','latex', ...
    'EdgeColor','none','FontSize',12,'FontWeight','bold');

% Optional: set a latex interpreter for the whole figure (some versions)
set(findall(fig,'-property','Interpreter'),'Interpreter','latex');

print(fig,'figure.png','-dpng','-r300')


