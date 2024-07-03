function [min_energy_consumed, min_energy_route , e_c ,tableData,ee_c,r, rr,ii,ccc,numRequests,x,minpower,status,table1] = calculateEnergyConsumption(n, cen1, s, d,dij,xij)
global tableData; %#ok<GVMIS>

l=0;
ivv=0;
Band=500;   
counter=0;
 Tlamda_sd=0;
lamdam=150;  
alpha=2;
consP=0;
Tlamda_sd=0;
   A_L_eq=[];   
    B_L_eq=[];    
    A_eq=[];
    B_eq=[];   

route_number = 0;
region_size1=1800;
 initial_node_energy = 500;
 node_energy = ones(1, n) * initial_node_energy;
  tableData = table('Size', [1, 5], 'VariableTypes', {'double', 'double', 'double', 'double', 'cell'}, 'VariableNames', {'request', 's', 'd', 'min_energy_consumed', 'min_energy_route'});


 while true
      numRequests=input(['Enter requests \n']);
for iter = 1:numRequests %#ok<ALIGN> 
    ivv=0;
     mido=[];
      r=[];
     m=0;

      l=l+1;
    tt=1;  
while tt ==1  
 disp('Enter the number of Source and Destination');
 disp('The value should not be more than the number of nodes');
 s=input('Source \n');
d=input('Destination \n');
if all(s<=n)&&all(d<=n)
    
tt=0;

end
% Generate the power vector
 pij=dij.^alpha; 
  PVij=pij';   
         PVij=PVij(:); 
        PVij=PVij';
        CC=PVij;  
  % Bandwidth constraint
% -----------------------------------

  AAA=[];   
for i=1:n   
    AA=zeros(n);  
    AA(i,:)=1;
    AA(:,i)=AA(:,i)+1;
    AA=AA';
    AA=AA(:);   
    AA=AA';
    AAA=[AAA;AA];   
end;    
    A_L_eq=[A_L_eq;AAA];  
    B_L_eq=[B_L_eq;Band/150*ones(n,1)];



   TPower=sum(sum(triu(pij).*triu(xij)));    
    A_L_eq=[A_L_eq;CC];    
    B_L_eq=[B_L_eq;TPower];   
      %  delta constraint
% -----------------------------------
       delta_sd=4;  
    A_L_eq=[A_L_eq;ones(1,n*n)];    
    B_L_eq=[B_L_eq;delta_sd]; 
%  routes constraint
%----------------------------------------
AA1=[];   
for i=1:n   
    AA=zeros(n);  
    AA(i,:)=1;
    AA(:,i)=AA(:,i)-1;
    AA=AA';
    AA=AA(:);   
    AA=AA';
    AA1=[AA1;AA];   
    if s==i
        BB=1;   
    elseif d==i   
        BB=-1;  
    else
        BB=0;   
    end;      
    B_eq=[B_eq;BB];         
end;    
A_eq=[A_eq;AA1];   
% Another set of less than or equal  
AA=eye(n*n);    
B=xij';
B=B(:);
A_L_eq=[A_L_eq;AA];    
B_L_eq=[B_L_eq;B]; 

 
       

%linear programming function
lb=zeros(1,n*n);   
ub=ones(1,n*n);   
e=2^-24;   
M=[1:n*n];
[x,minpower,status]=IP1(CC,A_L_eq,B_L_eq,A_eq,B_eq,lb,ub,M,e); 
if status ==1   
    xres=reshape(x,n,n);   
    xres=xres';   
    last=s;   
    pp1=[last]; 
    while (last ~= d)
        last=find(xres(last,:) == 1);
        pp1=[pp1,last];
    end;
end;    

if status ==1  
    path_string=[];  
    for k=1:length(pp1)-1
        path_string=[path_string,num2str(pp1(k)),'-->' ];
    end;
    path_string=[path_string,num2str(pp1(end))];    
    counter=counter+1;   
    consP=consP+minpower;  
    pij=pij+xres.*pij; 
%      Tlamda_sd=Tlamda_sd+lamda_sd; 

end
end
    m_lambda = 150; 
    requests = cell(n, 1);

        for i = 1:n
            k = poissrnd(1); 
            destinations = setdiff(1:n, i); 
            dest_indices = randsample(destinations, k, true);
            traffic_demands = normrnd(m_lambda, sqrt(0.5 * m_lambda), k, 1);

            requests{i} = [i * ones(k, 1), dest_indices', traffic_demands];
        end



    V = 1:n;
    OV = V;
    i_s = find(V ~= s);
    V = V(i_s);
    i_d = find(V ~= d);
    V = V(i_d);
    rr = [];
    c1 = 0;
% Delay constraint:
    for i = 1:length(V)
        cc = combnk(V, i);
        for ii = 1:size(cc, 1)
            pp = perms(cc(ii, :));
            for iii = 1:size(pp, 1)
                ccc = zeros(1, n - 2);
                ccc(1:size(pp, 2)) = pp(iii, :);
                rr = [rr; ccc];
            end
        end
    end

        route = [s d];
    if(dij(s,d) < region_size1) 
         m = m+ 1;
            r{m} = route;
    end
    
    for i = 1:length(rr)
        ivv=0;
          route=[];
        ii = find(rr(i, :) ~= 0);
        route_number = route_number + 1;
        route = [s rr(i, ii) d];
              x=[];
for ivv = 1:length(route)-1
            xxx=route(ivv);
            zzz=route(ivv + 1);
           x{ivv} = dij(xxx,zzz);
end

numeric_x = cell2mat(x);

if numeric_x <= region_size1
    m = m + 1;
    r{m} = route;
end    
    end
   


  %-----------------------------------------------------------\
   route_number = 0;
  
  for i=1:m 
   route_number = route_number + 1;
        energy_consumed = 0;
        nano=10^-6;
        Eamp =100*10^-9;
        b=1;
       mido=[r{i}];
     
         
        for j = 1:length(mido) - 1
            station1 = mido(j);
            station2 = mido(j + 1);
            distance = (sqrt(cen1(station1, 1) - cen1(station2, 1))^2 + (cen1(station1, 2) - cen1(station2, 2)));
              E_tx =Eamp*distance^2*b;
            node_energy(station1) = node_energy(station1) -E_tx;
             E_rx = 50 *nano * b; 
            node_energy(station2) = node_energy(station2) -E_rx;
            % Power constraint
            % The total power
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            energy_consumed=energy_consumed+((E_tx+E_rx)*10^-3);
            A_L_eq=[A_L_eq;CC];    
            B_L_eq=[B_L_eq;energy_consumed];   
        end
        eval(['route_number_', num2str(route_number), '=mido']); 
        disp(['Node energy levels after transmission: ', num2str(node_energy)]);
        disp(['Energy consumed for ', num2str(route_number-1), ': ', num2str(energy_consumed)]); 
  end 
    %---------------------------------------------
    
    min_energy_consumed = inf; % Initialize with a large value
    min_energy_route = [];
    e_c = zeros(1, route_number);
%Transmitting power constraint:
    for i = 1:route_number % Loop through each route
        route_var_name = ['route_number_', num2str(i)];
        current_route = eval(route_var_name);
             
        energy_consumed = 0;
         nano =10^-9;
        Eamp =100*10^-12;
        b=1;
        for j = 1:length(current_route) - 1
            station1 = current_route(j);
            station2 = current_route(j + 1);

            % Calculate distance between stations (you might use your dij matrix here)
            distance = (sqrt(cen1(station1, 1) - cen1(station2, 1))^2 + (cen1(station1, 2) - cen1(station2, 2)));

            % Calculate energy consumed based on the linear model
           E_tx =Eamp*distance^2*b;
           % transmation 
            node_energy(station1) = node_energy(station1) -E_tx;
             E_rx = 50 *nano*b; 
             %received
            node_energy(station2) = node_energy(station2) -E_rx;
            energy_consumed=energy_consumed+((E_tx+E_rx)*10^3);
          
        end
         e_c(i) = [energy_consumed]; 
%     [ee_c] = add_numbers(e_c);
            ee_c(iter ,:) ={e_c};
        if energy_consumed < min_energy_consumed
            min_energy_consumed = energy_consumed;
            min_energy_route = current_route;
        end
    end
    % Delay constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    % Display or store the total energy consumed for all routes
    disp(['the minimum energy consumption in all route ', num2str(min_energy_consumed)]);
    disp(['Route with Minimum Energy Consumption: ', num2str(min_energy_route)]);

    % Plot the energy consumption for each route
    figure;
    bar(e_c);
    xlabel('Route Number');
    ylabel('Energy Consumption (J)');
    title('Energy Consumption for Each Route');
    grid on;
saveas(gcf, fullfile('C:\Users\20101\OneDrive\Desktop\project lol\save_graphic', 'energy_consumption_bar_chart.jpg'));
    % Plot the energy consumption for each route
    figure;
    plot(1:1:route_number, e_c, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
    xlabel('Route Number');
    ylabel('Energy Consumption(J)');
    title('Energy Consumption for Each Route');
    grid on;
 saveas(gcf, fullfile('C:\Users\20101\OneDrive\Desktop\project lol\save_graphic', 'energy_consumption_line_chart.jpg'));
   
 
tableData(l, :) ={l, s, d, min_energy_consumed,{min_energy_route}};
disp(tableData); 
       
end 

  again = input('Do you need any Request again? yes/no (y/n)\n', 's');

    if strcmp(again, 'yes') || strcmp(again, 'y')
        
    else
        close(gcf); 
        break; 
    end
  end
    writetable(tableData, 'energy_consumption_data.csv');
        winopen('energy_consumption_data.csv');
  end
 
