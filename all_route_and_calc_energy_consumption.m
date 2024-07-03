clear all 
close all
clc
fclose('all');
try
    Excel = actxGetRunningServer('Excel.Application');
    Workbooks = Excel.Workbooks;
    Workbooks.Close;
    Excel.Quit;
    delete(Excel);
catch
    disp('No Excel instance found.');
end
format short g;
string_table=[];   
counter=0; 
tt=1;  
while tt ==1   
% ENTER NUMBER OF THE NOEDS 
n=input(['Enter the number of nodes you assumed,\n', ...
    'the value should be integer more than 1 \n']);  
n=round(n);
if n>1
tt=0;
end
end
%DISTRIBUTORED OF THE NODES
region_size1 = 1800;
xx = rand(1, n) * 3000;
yy = rand(1, n) * 3000;
  
for k=1:n
    xr=randi(length(xx));   
    xc(k)=xx(xr);  
    xi=find(xx ~= xc(k));
    xx=xx(xi);   
    yr=randi(length(yy));   
    yc(k)=yy(yr);  
    yi=find(yy ~= yc(k));
    yy=yy(yi);
end;    
%SET THE RANDOM OF THE NODES     
cen1=[xc',yc'];   
num1=1:n; 
stations_random_node(cen1,num1);  

dij=zeros(n,n);   
%CALCULATION OF THE DISTANCE BETWEEN I,J
for k=1:n   
    for kk=1:n   
        dij(k,kk)=sqrt((cen1(k,1)-cen1(kk,1))^2+(cen1(k,2)-cen1(kk,2))^2);
%          line([cen1(k, 1), cen1(kk, 1)], [cen1(k, 2), cen1(kk, 2)]);
    end;
end;
%DISCAVER THE NEIGHBOR AND ANGLE  
for k=1:n 
    ss=['station',num2str(k)];   
    TH=[];
    P1=cen1(k,:);
    neighbor=find( (dij(k,:) > 0) & (dij(k,:) <=region_size1));  
    if (length(neighbor) >=1)
        for kk=1:length(neighbor)
            P2=cen1(neighbor(kk),:);   
            theta=theta1(P1,P2);
            TH(kk)=theta;  
        end;
        eval([ss,'.neighborsNo = neighbor;']);   
        eval([ss,'.neighborsAng = TH;']);  
    else
        eval([ss,'.neighborsNo = [];']);   
        eval([ss,'.neighborsAng = [];']);  
    end;    
end;  
xij=zeros(n);
for k=1:n 
    ss=['station',num2str(k)];   
    k1=eval([ss,'.neighborsNo']);   
    xij(k,k1)=1;
end;       
m_lambda = 150; 
requests = cell(n, 1);

%Bandwidth constraint:
for i = 1:n
    k = poissrnd(1);
    destinations = setdiff(1:n, i);
    dest_indices = randsample(destinations, k, true);
    traffic_demands = normrnd(m_lambda, sqrt(0.5 * m_lambda), k, 1);
    requests{i} = [i * ones(k, 1), dest_indices', traffic_demands];
end
for i = 1:n
    for j = 1:size(requests{i}, 1)
        s = requests{i}(j, 1);
        d = requests{i}(j, 2);
%         plot([xc(s), xc(d)], [yc(s), yc(d)], 'b');
    end
end
scatter(xc, yc, 'filled');
hold on;
%Topology constraints:2&3
for k = 1:n
    ss = ['station', num2str(k)];
    neighbors = eval([ss, '.neighborsNo']);
    for kk = 1:length(neighbors)
        neighbor_index = neighbors(kk);
        plot([xc(k), xc(neighbor_index)], [yc(k), yc(neighbor_index)], 'r--');  
      
    end
end
hold off;
title('Nodes');
xlabel('X');
ylabel('Y');


[min_energy_consumed, min_energy_route, e_c, tableData,ee_c,r, rr,ii,ccc,numRequests,x,minpower,status] = calculateEnergyConsumption(n, cen1, s, d ,dij,xij);

num_routes = size(min_energy_route, 1);
num_items = size(e_c, 2);

% Repeat route numbers for each item
route_numbers = repmat((1:num_routes)', 1, num_items);  % Adjusted this line

% Reshape e_c into a column vector
energy_consumption_items = reshape(e_c', [], 1);

% Ensure both variables have the same number of rows
if numel(route_numbers) ~= numel(energy_consumption_items)
    error('Number of rows in route_numbers and energy_consumption_items must be the same.');
end

% Combine route number and energy consumption into a table  
i = 1:numel(energy_consumption_items);
route_table = table(i', energy_consumption_items, 'VariableNames', {'RouteNumber', 'EnergyConsumption'});  % Adjusted this line
  











 