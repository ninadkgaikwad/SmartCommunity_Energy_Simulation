function [U_k_PriorityStack] = PriorityStackController_SmartCommunity(E_LoadData,E_Mis)

% Author: Ninad Kiran Gaikwad
% Date: Mar/15/2021
% Description: PriorityStackController_SmartCommunity - PecanStreet Data PriorityStack Controller Logic

%% PriorityStackController_SmartCommunity - PecanStreet Data PriorityStack Controller Logic

%% Inital Basic Computation

% Getting Size of E_LoadData
[Row_E_LoadData,Column_E_LoadData, Depth_E_LoadData]=size(E_LoadData);

% Initializing U_k_PriorityStack
U_k_PriorityStack=zeros(1,Column_E_LoadData-9,Depth_E_LoadData);

% Getting absolute value of E_Mis
E_Mis_Abs=abs(E_Mis);

%% Priority Stack Controller Logic

% Getting Total Equipment wise energy usage for all houses
E_LoadData_Total=sum(E_LoadData(1,:,:),3);

% Getting E_Mis_Abs_Load_NotAC
E_Mis_Abs_Load_NotAC=E_Mis_Abs;

% Initializing E_Shedded
E_Shedded=0;

% Computing which Devices are shedded according to priority (Devices in E_LoadData are in Descending order of Priority)
for ii = Column_E_LoadData:-1:9+1 % For each device in Ascending order of Priority
    
    if (E_Shedded<E_Mis_Abs_Load_NotAC) % More Devices can be shedded
        
        % Updating E_Shedded
        E_Shedded = E_Shedded+E_LoadData_Total(ii);
        
        % Computing Load Control Commands (Load Shedded)
        U_k_PriorityStack(1,ii-9,1:Depth_E_LoadData)=zeros(1,1,Depth_E_LoadData);
        
        if (E_Shedded>=E_Mis_Abs_Load_NotAC)
            
            % Computing Load Control Commands (Load not Shedded)
            U_k_PriorityStack(1,1:(ii-9-1),1:Depth_E_LoadData)=ones(1,length(1:(ii-9-1)),Depth_E_LoadData);

            % Breaking from the For Loop
            break;
            
        end
  
    else % More Devices cannot be shedded 
        
        % Computing Load Control Commands (Load not Shedded)
        U_k_PriorityStack(1,1:(ii-9-1),1:Depth_E_LoadData)=ones(1,length(1:(ii-9-1)),Depth_E_LoadData);
        
        % Breaking from the For Loop
        break;
        
    end
    
end

   
end