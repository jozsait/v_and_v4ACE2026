clear all;
close all;
clc;

%% generate signals
N = 21;
t = linspace(0,2*pi,N);
f = sin(t);
noise = 0.2*randn(1,length(f));
f = f + noise;

%% reconstruction
order_array = 2:25;
error_array = order_array*0;
idx = 1;
for order=2:10
    t_r = linspace(0,2*pi,N);
    p_coeff = polyfit(t,f,order);
    f_r = polyval(p_coeff,t_r);
    
    disp({'the estimation error is ', mean(sqrt((f-f_r).^2)) })
    error_array(idx) = mean(sqrt((f-f_r).^2));
    idx = idx + 1;
end
%% visualisation
plot(t,f,'-r','LineWidth',2)
hold on
plot(t_r,f_r,'--k','LineWidth',2)

figure
plot(order_array,error_array,'ko',MarkerFaceColor='black')