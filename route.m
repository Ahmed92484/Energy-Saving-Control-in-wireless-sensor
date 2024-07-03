function [min_energy_consumed, min_energy_route] = route(s,d)
% s=input('Source \n');
% d=input('Destination \n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assume a value for the destination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d=4;

OV=V;
i_s=find(V~=s);
V=V(i_s);
i_d=find(V~=d);
V=V(i_d);
rr=[];   
c1=0;
for i=1:length(V)
    cc=combnk(V,i);
    for ii=1:size(cc,1)
        pp=perms(cc(ii,:));    
        for iii=1:size(pp,1)
            ccc=zeros(1,n-2);    
            ccc(1:size(pp,2))=pp(iii,:);
            rr=[rr;ccc];
        end; 
    end;
end;

% Remove the zeros and add the source and destination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

route_number=0;   
% The first route from source to destination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

route_number=route_number+1;   
route=[s d];
eval(['route_number_',num2str(route_number),'=route']);   
%disp(['Energy consumed for ',num2str(1), ': ', num2str(dij(s,d))]);

for i=1:size(rr,1)
    ii=find(rr(i,:) ~= 0);
    route_number=route_number+1;
    route=[s rr(i,ii) d]; 
    
    route_var_name = ['route_number_', num2str(i)];
    current_route = eval(route_var_name);
    energy_consumed=0;
    for j = 1:length(current_route)-1
        station1 = current_route(j);
        station2 = current_route(j+1);
        % Calculate distance between stations (you might use your dij matrix here)
        distance = sqrt((cen1(station1,1)-cen1(station2,1))^2 + (cen1(station1,2)-cen1(station2,2))^2);
        % Calculate energy consumed based on the linear model
        energy_consumed = energy_consumed +  distance ;
    end
    disp(['Energy consumed for ',num2str(route_number-1), ': ', num2str(energy_consumed)]);
    eval(['route_number_',num2str(route_number),'=route'])

end;
    for i = route_number:route_number% Loop through each route
        route_var_name = ['route_number_', num2str(i)];
        current_route = eval(route_var_name);
        
        energy_consumed = 0;
        for j = 1:length(current_route)-1
            station1 = current_route(j);
            station2 = current_route(j+1);
            
            % Calculate distance between stations (you might use your dij matrix here)
            distance = sqrt((cen1(station1,1)-cen1(station2,1))^2 + (cen1(station1,2)-cen1(station2,2))^2);
            
            % Calculate energy consumed based on the linear model
            energy_consumed = energy_consumed +  distance ;
        end;
        disp(['Energy consumed for ', route_var_name, ': ', num2str(energy_consumed)]);
    end
eval(['number_of_routes_for_n_eq_',num2str(n),'=route_number'])
%---------------------------------------------

min_energy_consumed = inf;  % Initialize with a large value
min_energy_route = [];
e_c=[];
for i = 1:route_number% Loop through each route
    route_var_name = ['route_number_', num2str(i)];
    current_route = eval(route_var_name);
    
    energy_consumed = 0;
    for j = 1:length(current_route)-1
        station1 = current_route(j);
        station2 = current_route(j+1);
        
        % Calculate distance between stations (you might use your dij matrix here)
        distance = sqrt((cen1(station1,1)-cen1(station2,1))^2 + (cen1(station1,2)-cen1(station2,2))^2);
        
        % Calculate energy consumed based on the linear model
        energy_consumed = energy_consumed +  distance ;
        e_c(i)=[energy_consumed];
    end
   if energy_consumed < min_energy_consumed
        min_energy_consumed = energy_consumed;
        min_energy_route = current_route;
    end
    
    % Display or store the results for each route
%     disp(['Energy consumed for ', route_var_name, ': ', num2str(energy_consumed)]);
end

% Display or store the total energy consumed for all routes
disp(['the minimum energy consumption in all route ', num2str(min_energy_consumed)]);
disp(['Route with Minimum Energy Consumption: ', num2str(min_energy_route)]);

% Plot the energy consumption for each route
figure;
bar(e_c);
xlabel('Route Number');
ylabel('Energy Consumption');
title('Energy Consumption for Each Route');
grid on;
% Display or store the total energy consumed for all routes
% disp(['the minimum energy consumption in all route ', num2str(min_energy_consumed)]);
% disp(['Route with Minimum Energy Consumption: ', num2str(min_energy_route)]);

% Plot the energy consumption for each route
figure;
plot(1:route_number, e_c, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Route Number');
ylabel('Energy Consumption');
title('Energy Consumption for Each Route');
grid on;


