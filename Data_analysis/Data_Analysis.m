%DATA ANALYSIS
clear
clc
close all
warning off
set(0,'defaultAxesFontSize',30)
set(0,'DefaultLineLineWidth',2)

%Raw data obtained from Star-CCM
%Description of the data can be found in Outputs_StarCCM.xlsx spreadsheet
Raw_Data = xlsread('Data.xlsx');

%Tabulating raw data
Table = table([NaN 0 0.2 0.4 0.8 1 -0.2 -0.4 -0.8 -1]',[30 45 60 90;Raw_Data(:,1:2:end)],[30 45 60 90;Raw_Data(:,2:2:end)]);

%Two main headers
%Columns:
%1) Actuation velocity, m/s
%2) D_bar, N (Drag - Drag_Actuation)
%3) Drag coefficient,C_d
%First row represents different angles of actuation with respect to horizontal, degrees
Table.Properties.VariableNames = {'Velocity','D_bar','C_d'};

for i=1:4
    Table.D_bar(2:end,i) = Table.D_bar(2,i)-Table.D_bar(2:end,i);
    Table.C_d(2:end,i) = Table.C_d(2,i)-Table.C_d(2:end,i);
end
disp('Data:');
disp(Table);

%Marker types
Marker = ['x' 'o' '*' 'p'];
%Linear fit colours
Colour = ['r' 'b' 'r' 'b'];
%Initialising loop
dD_dU_blow = zeros;
dCd_dU_blow = zeros;
for i = 1:4
    %Producing linear fit across data points
    fitobject = fit(Table.Velocity(2:5),Table.D_bar(2:5,i),'poly1');
    %Obtaining gradients
    coeffvals = coeffvalues(fitobject);
    dD_dU_blow(i) = coeffvals(1);
    
    %Plots
    %If statement splits up plots
    if i < 3
        figure(1)
    else
        figure(2)
    end
    plot(Table.Velocity(2:5),Table.D_bar(2:5,i),Marker(i),'Markersize',20,'Color',Colour(i))
    hold on
    plot(fitobject,Colour(i))
    grid on
    xlabel('$U_{act}$ (m/s)','Interpreter','Latex')
%     ylabel('D - D_{act}')
    ylabel('$\bar{D}$ (N)','Interpreter','Latex')
end

for i = 1:4
    %Producing linear fit across data points
    fitobject = fit(Table.Velocity(2:5),Table.C_d(2:5,i),'poly1');
    %Obtaining gradients
    coeffvals = coeffvalues(fitobject);
    dCd_dU_blow(i) = coeffvals(1);
    
    %Plots
    %If statement splits up plots
    if i < 3
        figure(3)
    else
        figure(4)
    end
    plot(Table.Velocity(2:5),Table.C_d(2:5,i),Marker(i),'Markersize',20,'Color',Colour(i))
    hold on
    plot(fitobject,Colour(i))
    grid on
    xlabel('$U_{act}$ (m/s)','Interpreter','Latex')
%     ylabel('C_d - C_{d_{act}}')
    ylabel('$\bar{C_d}$','Interpreter','Latex')
end

%Initialising loop
dD_dU_suck = zeros;
dCd_dU_suck = zeros;
for i = 1:4
    %Producing linear fit across data points
    fitobject = fit([Table.Velocity(2) Table.Velocity(7:end)']',[Table.D_bar(2,i) Table.D_bar(7:end,i)']','poly1');
    %Obtaining gradients
    coeffvals = coeffvalues(fitobject);
    dD_dU_suck(i) = coeffvals(1);
    
    %Plots
    %If statement splits up plots
    if i < 3
        figure(5)
    else
        figure(6)
    end
    plot([Table.Velocity(2) Table.Velocity(7:end)']',[Table.D_bar(2,i) Table.D_bar(7:end,i)']',Marker(i),'Markersize',20,'Color',Colour(i))
    hold on
    plot(fitobject,Colour(i))
    grid on
    xlabel('$U_{act}$ (m/s)','Interpreter','Latex')
%     ylabel('D - D_{act}')
    ylabel('$\bar{D}$ (N)','Interpreter','Latex')
end

for i = 1:4
    %Producing linear fit across data points
    fitobject = fit([Table.Velocity(2) Table.Velocity(7:end)']',[Table.C_d(2,i) Table.C_d(7:end,i)']','poly1');
    %Obtaining gradients
    coeffvals = coeffvalues(fitobject);
    dCd_dU_suck(i) = coeffvals(1);
    
    %Plots
    %If statement splits up plots
    if i < 3
        figure(7)
    else
        figure(8)
    end
    plot([Table.Velocity(2) Table.Velocity(7:end)']',[Table.C_d(2,i) Table.C_d(7:end,i)']',Marker(i),'Markersize',20,'Color',Colour(i))
    hold on
    plot(fitobject,Colour(i))
    grid on
    xlabel('$U_{act}$ (m/s)','Interpreter','Latex')
%     ylabel('C_d - C_{d_{act}}')
    ylabel('$\bar{C_d}$','Interpreter','Latex')
end

%Sensitivities
Output_blow = table(Table.D_bar(1,:)',dD_dU_blow',dCd_dU_blow');
Output_blow.Properties.VariableNames = {'Theta','dD_dU','dCd_dU'};
disp('Obtained sensitivities (BLOWING):');
disp(Output_blow);

Output_suck = table(Table.D_bar(1,:)',dD_dU_suck',dCd_dU_suck');
Output_suck.Properties.VariableNames = {'Theta','dD_dU','dCd_dU'};
disp('Obtained sensitivities (SUCKING):');
disp(Output_suck);

%Legends
figure(1)
legend('\theta = 30 degrees','{}','\theta = 45 degrees','{}','\theta = 60 degrees','{}','\theta = 90 degrees','')
figure(2)
legend('\theta = 60 degrees','{}','\theta = 90 degrees','')
figure(3)
legend('\theta = 30 degrees','{}','\theta = 45 degrees','{}','\theta = 60 degrees','{}','\theta = 90 degrees','')
figure(4)
legend('\theta = 60 degrees','{}','\theta = 90 degrees','')

figure(5)
legend('\theta = 30 degrees','{}','\theta = 45 degrees','{}','\theta = 60 degrees','{}','\theta = 90 degrees','')
figure(6)
legend('\theta = 60 degrees','{}','\theta = 90 degrees','')
figure(7)
legend('\theta = 30 degrees','{}','\theta = 45 degrees','{}','\theta = 60 degrees','{}','\theta = 90 degrees','')
figure(8)
legend('\theta = 60 degrees','{}','\theta = 90 degrees','')
