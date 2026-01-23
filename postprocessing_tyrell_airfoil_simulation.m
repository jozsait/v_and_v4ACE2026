clear all;
close all;
clc;


data_dir = 'C:\Users\Tamas.Jozsa\OneDrive - Cranfield University\Desktop\tyrell4ACE\postproc';
exp_data_file = 'Experiment_Cp_AoA_1_Freestream.dat';
sim_data_file = 'pressure_coeff_case_compr_turb_SA_Re4.60e5_farfield_BC_second_order_AUSM.dat';

exp_data = readmatrix(fullfile(data_dir, exp_data_file));
exp_data = array2table(exp_data, ...
    'VariableNames', {'x_coordinate', 'pressure_coefficient'});

sim_data = readtable(fullfile(data_dir, sim_data_file));
sim_data_upper = sortrows(sim_data(sim_data.pressure_coefficient >= 0, :), 'x_coordinate');
sim_data_lower = sortrows(sim_data(sim_data.pressure_coefficient <= 0, :), 'x_coordinate');

plot(exp_data.x_coordinate,exp_data.pressure_coefficient,'ko','MarkerFaceColor','black')
hold on
plot(sim_data_upper.x_coordinate,sim_data_upper.pressure_coefficient,'k-','LineWidth',2)
plot(sim_data_lower.x_coordinate,sim_data_lower.pressure_coefficient,'b--','LineWidth',2)
xlabel('x')
ylabel('C_p')
