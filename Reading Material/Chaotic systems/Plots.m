%Plots from Lorenz sys
clear 
clc
close all

load('Data.mat')

% 3D plot
figure
plot3(x,y,z);

%Plot of convergence of J
figure
plot(T_plot,J_mean_plot1)

%Plot of change in J with beta
figure
plot(beta_plot,J_mean_plot2(1,:),'r')
hold on
plot(beta_plot,J_mean_plot2(2,:),'b')
plot(beta_plot,J_mean_plot2(3,:),'g')