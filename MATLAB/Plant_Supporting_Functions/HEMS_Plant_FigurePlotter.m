function [ ] = HEMS_Plant_FigurePlotter(HEMS_Plant_FigurePlotter_Input)

% Author: Ninad Kiran Gaikwad
% Date: Mar/15/2021
% Description: HEMS_Plant_FigurePlotter_Input - Plotting Figures for the
% Baseline

close all;

EToP_Converter=(60/10); % kWh --> W

%% HEMS_Plant_FigurePlotter - Plotting Figures for Controller and Plant

%% Getting the desired data from the HEMS_Plant_Baseline_FigurePlotter_Input - Struct

%----------------HEMS_Plant_Baseline_FigurePlotter_Input------------------%
X_k_Plant_History=HEMS_Plant_FigurePlotter_Input.X_k_Plant_History;
U_k_History=HEMS_Plant_FigurePlotter_Input.U_k_History; 
E_LoadData=HEMS_Plant_FigurePlotter_Input.E_LoadData;
E_Load_Desired=HEMS_Plant_FigurePlotter_Input.E_Load_Desired;
HEMSWeatherData_Output=HEMS_Plant_FigurePlotter_Input.HEMSWeatherData_Output;
HEMSPlant_Params=HEMS_Plant_FigurePlotter_Input.HEMSPlant_Params;
Community_Params=HEMS_Plant_FigurePlotter_Input.Community_Params;
Baseline_Output_Images_Path=HEMS_Plant_FigurePlotter_Input.Baseline_Output_Images_Path;
Single_House_Plotting_Index=HEMS_Plant_FigurePlotter_Input.Single_House_Plotting_Index;
Simulation_Params=HEMS_Plant_FigurePlotter_Input.Simulation_Params;

%-----------------------HEMSWeatherData_Output----------------------------%
Ws=HEMSWeatherData_Output.Ws;
T_am=HEMSWeatherData_Output.T_am;
GHI=HEMSWeatherData_Output.GHI;
DNI=HEMSWeatherData_Output.DNI;
DateTimeVector=HEMSWeatherData_Output.DateTimeVector;
DateTime_Matrix=HEMSWeatherData_Output.DateTime_Matrix;

%---------------------------HEMSPlant_Params------------------------------%

E_AC=HEMSPlant_Params.E_AC;
T_AC_max=HEMSPlant_Params.T_AC_max;
T_AC_min=HEMSPlant_Params.T_AC_min;
ACLoad_StartUp_Power=HEMSPlant_Params.ACLoad_StartUp_Power;
Eff_Inv=HEMSPlant_Params.Eff_Inv;

Battery_Energy_Max=HEMSPlant_Params.Battery_Energy_Max;
Battery_Energy_Min=HEMSPlant_Params.Battery_Energy_Min;
MaxRate_Charging=HEMSPlant_Params.MaxRate_Charging;
MaxRate_Discharging=HEMSPlant_Params.MaxRate_Discharging;
Eff_Charging_Battery=HEMSPlant_Params.Eff_Charging_Battery;
Eff_Discharging_Battery=HEMSPlant_Params.Eff_Discharging_Battery;
MaxRate_Discharging_StartUp=HEMSPlant_Params.MaxRate_Discharging_StartUp;

%---------------------------Community_Params.-----------------------------%
N_House=Community_Params.N_House;
N_PV_Bat=Community_Params.N_PV_Bat;
N_Bat=Community_Params.N_Bat;
N_PV=Community_Params.N_PV;
N_None=Community_Params.N_None; 

%---------------------------Simulation_Params-----------------------------%

Simulation_StepSize=Simulation_Params.Simulation_StepSize;

%% Basic Computation and Generating plotable quantities

% House Numbers
N1=N_PV_Bat;
N2=N_Bat;
N3=N_PV;
N4=N_None;

% Truncating Plant History
X_k_Plant_History=X_k_Plant_History(1:end-1,:,:);

% PV Quantities - Individual Houses
House_PV_E_Available=X_k_Plant_History(:,1,:);
House_PV_E_Used=X_k_Plant_History(:,2,:);
House_PV_E_UnUsed=X_k_Plant_History(:,3,:);

% Battery Quantities - Individual Houses
House_Bat_E_State=X_k_Plant_History(:,4,:);
House_Bat_E_Charging=X_k_Plant_History(:,5,:);
House_Bat_E_Discharging=X_k_Plant_History(:,6,:);

% House Temperature Quantities - Individual Houses
House_Temprature=X_k_Plant_History(:,7,:);

% House Energy Quantities - Individual Houses
House_Bat_E_OtherLoad_Desired=E_Load_Desired(:,1,:);
House_Bat_E_ACLoad_Desired=U_k_History(:,3,:).*(E_AC/Eff_Inv);%***
House_Bat_E_TotalLoad_Desired=House_Bat_E_OtherLoad_Desired+House_Bat_E_ACLoad_Desired;

House_Bat_E_OtherLoad_Actual=X_k_Plant_History(:,12,:);
House_Bat_E_ACLoad_Actual=X_k_Plant_History(:,11,:);
House_Bat_E_TotalLoad_Actual=House_Bat_E_OtherLoad_Actual+House_Bat_E_ACLoad_Actual;

% House Controller Quantities - Individual Houses
House_Bat_Controller_Charging_Desired=U_k_History(:,1,:);
House_Bat_Controller_Discharging_Desired=U_k_History(:,2,:);

House_AC_Controller_TurnOn_Desired=U_k_History(:,3,:);
House_AC_Controller_TurnOn_Actual=X_k_Plant_History(:,21,:);


% PV Quantities - All Houses (Addup)
Community_PV_E_Available=sum(X_k_Plant_History(:,1,:),3);
Community_PV_E_Used=sum(X_k_Plant_History(:,2,:),3);
Community_PV_E_UnUsed=sum(X_k_Plant_History(:,3,:),3);

% Battery Quantities - Individual Houses - All Houses (Addup)
Community_Bat_E_State=sum(X_k_Plant_History(:,4,:),3);
Community_Bat_E_Charging=sum(X_k_Plant_History(:,5,:),3);
Community_Bat_E_Discharging=sum(X_k_Plant_History(:,6,:),3);

% House Temperature Quantities - Individual Houses - All Houses (Average)
Community_Temprature=mean(X_k_Plant_History(:,7,:),3);

% House Energy Quantities - Individual Houses - All Houses (Addup)
Community_Bat_E_OtherLoad_Desired=sum(E_Load_Desired(:,1,:),3);
Community_Bat_E_ACLoad_Desired=sum(House_Bat_E_ACLoad_Desired(:,1,:),3);
Community_Bat_E_TotalLoad_Desired=Community_Bat_E_ACLoad_Desired+Community_Bat_E_OtherLoad_Desired;

Community_Bat_E_OtherLoad_Actual=sum(X_k_Plant_History(:,12,:),3);
Community_Bat_E_ACLoad_Actual=sum(X_k_Plant_History(:,11,:),3);
Community_Bat_E_TotalLoad_Actual=Community_Bat_E_OtherLoad_Actual+Community_Bat_E_ACLoad_Actual;

% House Controller Quantities - Individual Houses - All Houses (Addup)
Community_Bat_Controller_Charging_Desired=sum(U_k_History(:,1,:),3);
Community_Bat_Controller_Discharging_Desired=sum(U_k_History(:,2,:),3);

Community_AC_Controller_TurnOn_Desired=sum(U_k_History(:,3,:),3);
Community_AC_Controller_TurnOn_Actual=sum(X_k_Plant_History(:,21,:),3);

% Computing House and Community Battery SoC
House_Bat_SoC=(((House_Bat_E_State)-Battery_Energy_Min)/(Battery_Energy_Max-Battery_Energy_Min))*100;
Community_Bat_SoC=(((Community_Bat_E_State)-((N1+N2)*Battery_Energy_Min))/(((N1+N2)*Battery_Energy_Max)-((N1+N2)*Battery_Energy_Min)))*100;

% Computing House and Community Generation and Demand
House_E_Generation=House_PV_E_Used+House_Bat_E_Discharging;
House_E_Demand=House_Bat_E_Charging+House_Bat_E_TotalLoad_Actual;

Community_E_Generation=Community_PV_E_Used+Community_Bat_E_Discharging;
Community_E_Demand=Community_Bat_E_Charging+Community_Bat_E_TotalLoad_Actual;

% House/Community Level Battery Charging_Dispatchable/Discharging_Dispatchable
House_Bat_E_Charging_Dispatchable=zeros(length(GHI),1,N_House); % Initialization
House_Bat_E_Discharging_Dispatchable=zeros(length(GHI),1,N_House); % Initialization
for ii=1:length(GHI)
    if (Community_PV_E_Available(ii,1,1)>=Community_Bat_E_TotalLoad_Desired(ii,1,1))

        for jj=1:(N_PV_Bat+N_Bat)
            
            House_Bat_E_Charging_Dispatchable(ii,1,jj)=min(MaxRate_Charging*Simulation_StepSize,(Battery_Energy_Max-(X_k_Plant_History(ii,4,jj)))/Eff_Charging_Battery);
            House_Bat_E_Discharging_Dispatchable(ii,1,jj)=0;
            
        end

    else

        for jj=1:(N_PV_Bat+N_Bat)
            
            House_Bat_E_Charging_Dispatchable(ii,1,jj)=0;
            House_Bat_E_Discharging_Dispatchable(ii,1,jj)=min(MaxRate_Discharging*Simulation_StepSize,(X_k_Plant_History(ii,4,jj)-Battery_Energy_Min)*Eff_Discharging_Battery);
            
        end

    end
    
end

Community_Bat_E_Charging_Dispatchable=sum(House_Bat_E_Charging_Dispatchable(:,1,:),3);%***
Community_Bat_E_Discharging_Dispatchable=sum(House_Bat_E_Discharging_Dispatchable(:,1,:),3);%***

% House AC Startup Power Quantities
Community_AC_P_StartUp_Available=(Community_PV_E_Available/Simulation_StepSize)+(Community_Bat_Controller_Charging_Desired*MaxRate_Discharging_StartUp);

House_AC_StartUp_Desired=House_AC_Controller_TurnOn_Desired(:,1,:)-X_k_Plant_History(:,30,:);
House_AC_StartUp_Desired = (abs(House_AC_StartUp_Desired)+House_AC_StartUp_Desired)/2;
Community_AC_StartUp_Desired = sum(House_AC_StartUp_Desired(:,1,:),3);
Community_AC_P_StartUp_Required=Community_AC_StartUp_Desired*ACLoad_StartUp_Power;

House_AC_StartUp_Actual=X_k_Plant_History(:,21,:)-X_k_Plant_History(:,30,:);
House_AC_StartUp_Actual = (abs(House_AC_StartUp_Actual)+House_AC_StartUp_Actual)/2;
Community_AC_StartUp_Actual = sum(House_AC_StartUp_Actual(:,1,:),3);
Community_AC_P_StartUp_Used=Community_AC_StartUp_Actual*ACLoad_StartUp_Power;

% Creating Grouping Indices
N_All_Indices=1:N_House;
N_PV_Bat_Only_Indices=1:N1;
N_Bat_Only_Indices=N1+1:N2;
N_PV_Only_Indices=N2+1:N3;
N_None_Only_Indices=N3+1:N4;

%% DateTime Vector to Hours

D=DateTimeVector(2)-DateTimeVector(1);

M=minutes(D);

H=M/60;

L=length(DateTimeVector);

HoursVector=zeros(L,1);
HoursVector=zeros(L,1);

for ii=2:L
    HoursVector(ii,1)=HoursVector(ii-1,1)+H;
end

Len_Hours_Vector=length(HoursVector);

%% Plotting the Figures - Individual Houses and Community

% House PV Power Plots - Availabe/Used/Unused
h1=figure(1);
hold on
box on

for jj=union(N_PV_Bat_Only_Indices,N_PV_Only_Indices)
    
    P2 = plot(HoursVector,EToP_Converter*House_PV_E_Available(1:Len_Hours_Vector,1,jj),'-k','LineWidth',0.3);
    P3 = plot(HoursVector,EToP_Converter*House_PV_E_Used(1:Len_Hours_Vector,1,jj),'-b','LineWidth',0.4);
    P4 = plot(HoursVector,EToP_Converter*House_PV_E_UnUsed(1:Len_Hours_Vector,1,jj),'-r','LineWidth',0.2);

end
    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('House Level - PV Power','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3,P4],'PV Power Available','PV Power Used','PV Power Unused');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h1,strcat(Baseline_Output_Images_Path,'House_PV_Power_A_U_UnU'),'jpeg');

% Community PV Power Plots - Availabe/Used/Unused
h2=figure(2);
hold on
box on

P2 = plot(HoursVector,EToP_Converter*Community_PV_E_Available(1:Len_Hours_Vector,1,1),'-k','LineWidth',1.5);
P3 = plot(HoursVector,EToP_Converter*Community_PV_E_Used(1:Len_Hours_Vector,1,1),'-b','LineWidth',2);
P4 = plot(HoursVector,EToP_Converter*Community_PV_E_UnUsed(1:Len_Hours_Vector,1,1),'-r','LineWidth',1);
    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('Community Level - PV Power','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3,P4],'PV Power Available','PV Power Used','PV Power Unused');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h2,strcat(Baseline_Output_Images_Path,'Community_PV_Power_A_U_UnU'),'jpeg');

% House Battery SoC/Bat_Charging/Bat_Discharging Plots 
h3=figure(3);
hold on
box on

yyaxis left

for jj=union(N_PV_Bat_Only_Indices,N_Bat_Only_Indices)
    
    P2 = plot(HoursVector,House_Bat_SoC(1:Len_Hours_Vector,1,jj),'-b','LineWidth',0.3);

end
P3 = plot(HoursVector,100*ones(Len_Hours_Vector,1),'--k','LineWidth',1);
P4 = plot(HoursVector,0*ones(Len_Hours_Vector,1),'--k','LineWidth',1);

    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-5 115]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('$\%$','Interpreter','latex','FontSize', 14);
title('House Level - Battery SoC/Charging/Discharging','Interpreter','latex','FontSize', 14);

yyaxis right

for jj=union(N_PV_Bat_Only_Indices,N_Bat_Only_Indices)
    
    P5 = plot(HoursVector,EToP_Converter*House_Bat_E_Charging(1:Len_Hours_Vector,1,jj),'-g','LineWidth',0.3);
    P6 = plot(HoursVector,EToP_Converter*House_Bat_E_Discharging(1:Len_Hours_Vector,1,jj),'-r','LineWidth',0.3);

end
    
%ylim([-5 115]);
xlim([0 170]);
ax1 = ancestor(P4, 'axes');
yrule_2 = ax1.YAxis;
yrule_2(2).FontSize=14; 
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P5,P6],'SoC','Battery Charging Power','Battery Discharging Power');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h3,strcat(Baseline_Output_Images_Path,'House_Bat_SoC_C_DisC'),'jpeg');


% Community Battery SoC/Bat_Charging/Bat_Discharging Plots 
h4=figure(4);
hold on
box on



P2 = plot(HoursVector,Community_Bat_SoC(1:Len_Hours_Vector,1,1),'-b','LineWidth',1.5);
P3 = plot(HoursVector,100*ones(Len_Hours_Vector,1),'--k','LineWidth',1);
P4 = plot(HoursVector,0*ones(Len_Hours_Vector,1),'--k','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-5 115]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('SoC','Interpreter','latex','FontSize', 14);
title('Community Level - Battery SoC/Charging/Discharging','Interpreter','latex','FontSize', 14);

yyaxis right
   
P5 = plot(HoursVector,EToP_Converter*Community_Bat_E_Charging(1:Len_Hours_Vector,1,1),'-g','LineWidth',1);
P6 = plot(HoursVector,EToP_Converter*Community_Bat_E_Discharging(1:Len_Hours_Vector,1,1),'-r','LineWidth',1);

%ylim([-5 115]);
xlim([0 170]);
ax1 = ancestor(P4, 'axes');
yrule_2 = ax1.YAxis;
yrule_2(2).FontSize=14; 
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P5,P6],'SoC','Battery Charging Power','Battery Discharging Power');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h4,strcat(Baseline_Output_Images_Path,'Community_Bat_SoC_C_DisC'),'jpeg');

% House Temperature
h5=figure(5);
hold on
box on

for jj=1:N_House
    
    P2 = plot(HoursVector,House_Temprature(1:Len_Hours_Vector,1,jj),'-b','LineWidth',0.2);

end
P3 = plot(HoursVector,T_AC_max*ones(Len_Hours_Vector,1),'--k','LineWidth',1);
P4 = plot(HoursVector,T_AC_min*ones(Len_Hours_Vector,1),'--k','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Temperature $(^{\circ}C)$','Interpreter','latex','FontSize', 14);
title('House Level - Temperature','Interpreter','latex','FontSize', 14);

%legend1=legend([P2],'PV Power Available','PV Power Used','PV Power Unused');
%set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h5,strcat(Baseline_Output_Images_Path,'House_Temperature'),'jpeg');

% Community Temperature
h6=figure(6);
hold on
box on
    
P2 = plot(HoursVector,Community_Temprature(1:Len_Hours_Vector,1,1),'-b','LineWidth',1.5);
P3 = plot(HoursVector,T_AC_max*ones(Len_Hours_Vector,1),'--k','LineWidth',1);
P4 = plot(HoursVector,T_AC_min*ones(Len_Hours_Vector,1),'--k','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Temperature $(^{\circ}C)$','Interpreter','latex','FontSize', 14);
title('Community Level - Temperature','Interpreter','latex','FontSize', 14);

%legend1=legend([P2],'PV Power Available','PV Power Used','PV Power Unused');
%set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h6,strcat(Baseline_Output_Images_Path,'Community_Temperature'),'jpeg');

% House Load Power - Load/AC/Total Desired/Actual
h7=figure(7);

subplot(3,1,1)
hold on
box on

for jj=1:N_House
    
    P2 = plot(HoursVector,EToP_Converter*House_Bat_E_OtherLoad_Desired(1:Len_Hours_Vector,1,jj),'--b','LineWidth',0.3);
    P5 = plot(HoursVector,EToP_Converter*House_Bat_E_OtherLoad_Actual(1:Len_Hours_Vector,1,jj),'-b','LineWidth',0.2);

end

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('House Level - Load Power Desired/Actual','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P5],'Other Desired','Other Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;
  
subplot(3,1,2)
hold on
box on

for jj=1:N_House
    
    P3 = plot(HoursVector,EToP_Converter*House_Bat_E_ACLoad_Desired(1:Len_Hours_Vector,1,jj),'--r','LineWidth',0.3);
    P6 = plot(HoursVector,EToP_Converter*House_Bat_E_ACLoad_Actual(1:Len_Hours_Vector,1,jj),'-r','LineWidth',0.2);

end

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P3, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P3,P6],'AC Desired','AC Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

subplot(3,1,3)
hold on
box on

for jj=1:N_House
    
    P4 = plot(HoursVector,EToP_Converter*House_Bat_E_TotalLoad_Desired(1:Len_Hours_Vector,1,jj),'--k','LineWidth',0.3);
    P7 = plot(HoursVector,EToP_Converter*House_Bat_E_TotalLoad_Actual(1:Len_Hours_Vector,1,jj),'-k','LineWidth',0.2);

end

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P4, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P4,P7],'Total Desired','Total Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h7,strcat(Baseline_Output_Images_Path,'House_Loads_Other_AC_Total'),'jpeg');

% Community Load Power - Load/AC/Total Desired/Actual
h8=figure(8);

subplot(3,1,1)
hold on
box on

P2 = plot(HoursVector,EToP_Converter*Community_Bat_E_OtherLoad_Desired(1:Len_Hours_Vector,1,1),'--b','LineWidth',1.5);
P5 = plot(HoursVector,EToP_Converter*Community_Bat_E_OtherLoad_Actual(1:Len_Hours_Vector,1,1),'-b','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('Community Level - Load Power Desired/Actual','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P5],'Other Desired','Other Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

subplot(3,1,2)
hold on
box on

P3 = plot(HoursVector,EToP_Converter*Community_Bat_E_ACLoad_Desired(1:Len_Hours_Vector,1,1),'--r','LineWidth',1.5);
P6 = plot(HoursVector,EToP_Converter*Community_Bat_E_ACLoad_Actual(1:Len_Hours_Vector,1,1),'-r','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P3, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P3,P6],'AC Desired','AC Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

subplot(3,1,3)
hold on
box on

P4 = plot(HoursVector,EToP_Converter*Community_Bat_E_TotalLoad_Desired(1:Len_Hours_Vector,1,1),'--k','LineWidth',1.5);
P7 = plot(HoursVector,EToP_Converter*Community_Bat_E_TotalLoad_Actual(1:Len_Hours_Vector,1,1),'-k','LineWidth',1);
    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P4, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P4,P7],'Total Desired','Total Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h8,strcat(Baseline_Output_Images_Path,'Community_Loads_Other_AC_Total'),'jpeg');

% House Battery Controller - Charging/Discharging
h9=figure(9);
hold on
box on


for jj=union(N_PV_Bat_Only_Indices,N_Bat_Only_Indices)
    
    P2 = plot(HoursVector,House_Bat_Controller_Charging_Desired(1:Len_Hours_Vector,1,jj),'-r','LineWidth',0.3);
    P3 = plot(HoursVector,House_Bat_Controller_Discharging_Desired(1:Len_Hours_Vector,1,jj),'--b','LineWidth',0.3);
    
end
    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
ylim([0 2]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Battery Commands','Interpreter','latex','FontSize', 14);
title('House Level - Battery Controller Commands','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3],'Charging Command','Discharging Command');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h9,strcat(Baseline_Output_Images_Path,'House_Bat_Controller_Charging_Discharging'),'jpeg');

% Community Battery Controller - Charging/Discharging
h10=figure(10);
hold on
box on

P2 = plot(HoursVector,Community_Bat_Controller_Charging_Desired(1:Len_Hours_Vector,1,1),'-r','LineWidth',1);
P3 = plot(HoursVector,Community_Bat_Controller_Discharging_Desired(1:Len_Hours_Vector,1,1),'--b','LineWidth',1);
  
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([0 2]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Battery Commands','Interpreter','latex','FontSize', 14);
title('Community Level - Battery Controller Commands','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3],'Charging Command','Discharging Command');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h10,strcat(Baseline_Output_Images_Path,'Community_Bat_Controller_Charging_Discharging'),'jpeg');


% House AC Controller - on-off Desired/on-off Actual
h11=figure(11);
hold on
box on


for jj=1:N_House
    
    P2 = plot(HoursVector,House_AC_Controller_TurnOn_Desired(1:Len_Hours_Vector,1,jj),'-r','LineWidth',0.2);
    P3 = plot(HoursVector,House_AC_Controller_TurnOn_Actual(1:Len_Hours_Vector,1,jj),'--b','LineWidth',0.2);
    
end
    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
ylim([0 2]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('AC Commands','Interpreter','latex','FontSize', 14);
title('House Level - AC Controller Commands','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3],'Desired Command','Actual Command');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h11,strcat(Baseline_Output_Images_Path,'House_AC_Controller_Desired_Actual'),'jpeg');

% Community AC Controller - on-off Desired/on-off Actual
h12=figure(12);
hold on
box on

P2 = plot(HoursVector,Community_AC_Controller_TurnOn_Desired(1:Len_Hours_Vector,1,1),'-r','LineWidth',1);
P3 = plot(HoursVector,Community_AC_Controller_TurnOn_Actual(1:Len_Hours_Vector,1,1),'--b','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([0 2]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('AC Commands','Interpreter','latex','FontSize', 14);
title('Community Level - AC Controller Commands','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3],'Desired Command','Actual Command');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h12,strcat(Baseline_Output_Images_Path,'Community_AC_Controller_Desired_Actual'),'jpeg');

% Community AC Startup Power
h13=figure(13);
hold on
box on

P2 = plot(HoursVector,Community_AC_P_StartUp_Available(1:Len_Hours_Vector,1,1),'-k','LineWidth',1.5);
P3 = plot(HoursVector,Community_AC_P_StartUp_Required(1:Len_Hours_Vector,1,1),'-r','LineWidth',1);
P4 = plot(HoursVector,Community_AC_P_StartUp_Used(1:Len_Hours_Vector,1,1),'-b','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([0 2]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('House Level - AC Startup Power','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3,P4],'Start-Up Available','Start-Up Required','Start-Up Used');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h13,strcat(Baseline_Output_Images_Path,'Community_AC_StartUpPower_Av_Req_U'),'jpeg');

%% Plotting the Figures - Single House

Single_House_Plotting_Index=1;

% House PV Power Plots - Availabe/Used/Unused
h14=figure(14);
hold on
box on

P2 = plot(HoursVector,EToP_Converter*House_PV_E_Available(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-k','LineWidth',1.5);
P3 = plot(HoursVector,EToP_Converter*House_PV_E_Used(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-b','LineWidth',2);
P4 = plot(HoursVector,EToP_Converter*House_PV_E_UnUsed(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-r','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('House Level - PV Power','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3,P4],'PV Power Available','PV Power Used','PV Power Unused');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h14,strcat(Baseline_Output_Images_Path,'SingleHouse_PV_Power_A_U_UnU'),'jpeg');

% House Battery SoC/Bat_Charging/Bat_Discharging Plots 
h15=figure(15);
hold on
box on

yyaxis left

P2 = plot(HoursVector,House_Bat_SoC(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-b','LineWidth',1.5);

P3 = plot(HoursVector,100*ones(Len_Hours_Vector,1),'--k','LineWidth',1);
P4 = plot(HoursVector,0*ones(Len_Hours_Vector,1),'--k','LineWidth',1);

    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-5 115]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('$\%$','Interpreter','latex','FontSize', 14);
title('House Level - Battery SoC/Charging/Discharging','Interpreter','latex','FontSize', 14);

yyaxis right

P5 = plot(HoursVector,EToP_Converter*House_Bat_E_Charging(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-g','LineWidth',1);
P6 = plot(HoursVector,EToP_Converter*House_Bat_E_Discharging(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-r','LineWidth',1);

%ylim([-5 115]);
xlim([0 170]);
ax1 = ancestor(P4, 'axes');
yrule_2 = ax1.YAxis;
yrule_2(2).FontSize=14; 
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P5,P6],'SoC','Battery Charging Power','Battery Discharging Power');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h15,strcat(Baseline_Output_Images_Path,'SingleHouse_Bat_SoC_C_DisC'),'jpeg');

% House Temperature
h16=figure(16);
hold on
box on

P2 = plot(HoursVector,House_Temprature(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-b','LineWidth',1.5);

P3 = plot(HoursVector,T_AC_max*ones(Len_Hours_Vector,1),'--k','LineWidth',1);
P4 = plot(HoursVector,T_AC_min*ones(Len_Hours_Vector,1),'--k','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Temperature $(^{\circ}C)$','Interpreter','latex','FontSize', 14);
title('House Level - Temperature','Interpreter','latex','FontSize', 14);

%legend1=legend([P2],'PV Power Available','PV Power Used','PV Power Unused');
%set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h16,strcat(Baseline_Output_Images_Path,'SingleHouse_Temperature'),'jpeg');

% House Load Power - Load/AC/Total Desired/Actual
h17=figure(17);

subplot(3,1,1)
hold on
box on

P2 = plot(HoursVector,EToP_Converter*House_Bat_E_OtherLoad_Desired(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'--b','LineWidth',1);
P5 = plot(HoursVector,EToP_Converter*House_Bat_E_OtherLoad_Actual(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-b','LineWidth',1.5);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('House Level - Load Power Desired/Actual','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P5],'Other Desired','Other Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

subplot(3,1,2)
hold on
box on

P3 = plot(HoursVector,EToP_Converter*House_Bat_E_ACLoad_Desired(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'--r','LineWidth',1);
P6 = plot(HoursVector,EToP_Converter*House_Bat_E_ACLoad_Actual(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-r','LineWidth',1.5);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P3, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P3,P6],'AC Desired','AC Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

subplot(3,1,3)
hold on
box on

P4 = plot(HoursVector,EToP_Converter*House_Bat_E_TotalLoad_Desired(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'--k','LineWidth',1);
P7 = plot(HoursVector,EToP_Converter*House_Bat_E_TotalLoad_Actual(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-k','LineWidth',1.5);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P4, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

legend1=legend([P4,P7],'Total Desired','Total Actual');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h17,strcat(Baseline_Output_Images_Path,'SingleHouse_Loads_Other_AC_Total'),'jpeg');

% House Battery Controller - Charging/Discharging
h18=figure(18);
hold on
box on

P2 = plot(HoursVector,House_Bat_Controller_Charging_Desired(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'-r','LineWidth',1);
P3 = plot(HoursVector,House_Bat_Controller_Discharging_Desired(1:Len_Hours_Vector,1,Single_House_Plotting_Index),'--b','LineWidth',1);

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
ylim([0 2]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Battery Commands','Interpreter','latex','FontSize', 14);
title('House Level - Battery Controller Commands','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3],'Charging Command','Discharging Command');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h18,strcat(Baseline_Output_Images_Path,'SingleHouse_Bat_Controller_Charging_Discharging'),'jpeg');

% House AC Controller - on-off Desired/on-off Actual
h19=figure(19);
hold on
box on

for jj=1:N_House
    
    P2 = plot(HoursVector,House_AC_Controller_TurnOn_Desired(1:Len_Hours_Vector,1,jj),'-r','LineWidth',0.2);
    P3 = plot(HoursVector,House_AC_Controller_TurnOn_Actual(1:Len_Hours_Vector,1,jj),'--b','LineWidth',0.2);
    
end
    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
ylim([0 2]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('AC Commands','Interpreter','latex','FontSize', 14);
title('House Level - AC Controller Commands','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3],'Desired Command','Actual Command');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h19,strcat(Baseline_Output_Images_Path,'SingleHouse_AC_Controller_Desired_Actual'),'jpeg');


%% %% Plotting the Figures - Individual Houses and Community

House_E_Generation=House_PV_E_Used+House_Bat_E_Discharging;
House_E_Demand=House_Bat_E_Charging+House_Bat_E_TotalLoad_Actual;

Community_E_Generation=Community_PV_E_Used+Community_Bat_E_Discharging;
Community_E_Demand=Community_Bat_E_Charging+Community_Bat_E_TotalLoad_Actual;


% House Power Plots - Generation/Demand
h20=figure(20);
hold on
box on

for jj=1:N_House
    
    P2 = plot(HoursVector,EToP_Converter*House_E_Generation(1:Len_Hours_Vector,1,jj),'-b','LineWidth',0.3);
    P3 = plot(HoursVector,EToP_Converter*House_E_Demand(1:Len_Hours_Vector,1,jj),'--r','LineWidth',0.2);
    
    P4 = plot(HoursVector,EToP_Converter*(House_PV_E_Available(1:Len_Hours_Vector,1,jj)+House_Bat_E_Discharging_Dispatchable(1:Len_Hours_Vector,1,jj)),'-.b','LineWidth',0.3);
    P5 = plot(HoursVector,EToP_Converter*(House_Bat_E_TotalLoad_Desired(1:Len_Hours_Vector,1,jj)+House_Bat_E_Charging_Dispatchable(1:Len_Hours_Vector,1,jj)),'-.r','LineWidth',0.2);    

end


    
xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('House Level - Generation/Demand','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3,P4,P5],'Generation','Demand','Available Generation','Demand Desired');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h20,strcat(Baseline_Output_Images_Path,'House_P_Generated_Demand'),'jpeg');

% Community Power Plots - Generation/Demand
h21=figure(21);
hold on
box on

P2 = plot(HoursVector,EToP_Converter*Community_E_Generation(1:Len_Hours_Vector,1,1),'-b','LineWidth',1.5);
P3 = plot(HoursVector,EToP_Converter*Community_E_Demand(1:Len_Hours_Vector,1,1),'--r','LineWidth',1);

P4 = plot(HoursVector,EToP_Converter*(Community_PV_E_Available(1:Len_Hours_Vector,1,1)+Community_Bat_E_Discharging_Dispatchable(1:Len_Hours_Vector,1,1)),'-.b','LineWidth',0.3);
P5 = plot(HoursVector,EToP_Converter*(Community_Bat_E_TotalLoad_Desired(1:Len_Hours_Vector,1,1)+Community_Bat_E_Charging_Dispatchable(1:Len_Hours_Vector,1,1)),'-.r','LineWidth',0.2);    

xticks([0 24 48 72 96 120 144 168])
xticklabels({'0','24','48','72','96','120','144','168'})
%ylim([-6 25]);
xlim([0 170]);
xlabel('Time ($hours$)','Interpreter','latex','FontSize', 14);
ax1 = ancestor(P2, 'axes');
xrule_1 = ax1.XAxis;
xrule_1.FontSize=14;
yrule_1 = ax1.YAxis;
yrule_1(1).FontSize=14;
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);
title('Community Level - Generation/Demand','Interpreter','latex','FontSize', 14);

legend1=legend([P2,P3,P4,P5],'Generation','Demand','Available Generation','Demand Desired');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h21,strcat(Baseline_Output_Images_Path,'Community_P_Generated_Demand'),'jpeg');

%%

% House Battery SoC/Bat_Charging/Bat_Discharging Plots 
h22=figure(22);
hold on
box on

yyaxis right

for jj=union(N_PV_Bat_Only_Indices,N_Bat_Only_Indices)
    
    P5 = plot(HoursVector,EToP_Converter*House_Bat_E_Charging(1:Len_Hours_Vector,1,jj),'-g','LineWidth',0.3);
    P6 = plot(HoursVector,EToP_Converter*House_Bat_E_Discharging(1:Len_Hours_Vector,1,jj),'-r','LineWidth',0.3);

end
    
%ylim([-5 115]);
xlim([0 170]);
ax1 = ancestor(P5, 'axes');
yrule_2 = ax1.YAxis;
yrule_2(2).FontSize=14; 
ylabel('Power $(kW)$','Interpreter','latex','FontSize', 14);

title('House Level - Charging','Interpreter','latex','FontSize', 14);

%legend1=legend([P5],'SoC','Battery Charging Power','Battery Discharging Power');
set(legend1,'Interpreter','latex','FontSize', 12);

box off
hold off;

saveas(h22,strcat(Baseline_Output_Images_Path,'House_Bat_C'),'jpeg');



end

