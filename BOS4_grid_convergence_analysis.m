clear all;
close all;
clc;

%% simulation results
% domain size [m^2]
A = 76;

% number of cells
N = [inf, 18000, 8000, 4500];

% reattachment length
L = [0, 6.063, 5.972, 5.863];

% order of convergence
p = 1.53;

%% grid convergence computations
% representative grid spacing
h = sqrt(A./N);

% refinement ratios
r32 = h(4)/h(3);
r21 = h(3)/h(2);

% Richardson extrapolation
L(1) = (r21^p*L(2)-L(3))/(r21^p-1);

% relative change
eps21 = abs((L(3)-L(2))/L(2));
eps32 = abs((L(4)-L(3))/L(3));

% grid convergence indices
GCI21 = 1.25*eps21/(r21^p-1);
GCI32 = 1.25*eps32/(r32^p-1);

AR = r21^p*GCI21/GCI32;

%% visualisation
plot(h(2:end)/h(2),L(2:end),'-ko','MarkerFaceColor','black',LineWidth=2)
hold on
plot(h(1)/h(2),L(1),'rd','MarkerFaceColor','red','MarkerSize',8)
plot(h(1:2)/h(2),L(1:2),'--k',LineWidth=2)

xlabel('Normalised grid spacing')
ylabel('Rettachment length')