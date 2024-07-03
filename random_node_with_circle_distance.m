clear all
close all
clc
% Assume minimum distance in X, and Y between
% stations
string_table=[];   
% Assume number of stations n
tt=1;  
while tt ==1   
         n=input(['Enter the number of nodes you assumed,\n', ...
             'the value should be integer more than 1 \n']);  

         n=round(n);
     if n>1
           tt=0;
     end
end

region_size1 = 1800;
xc = rand(1, n) * region_size1;
yc = rand(1, n) * region_size1;
  
% Random locations selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
cen1=[xc',yc'];   
num1=1:n; 
stations_random_node(cen1,num1);  

dij=zeros(n,n);   

for k=1:n   
    for kk=1:n   
        dij(k,kk)=sqrt((cen1(k,1)-cen1(kk,1))^2+(cen1(k,2)-cen1(kk,2))^2);    
            % if dij(k, kk) <= 700
                
             %  line([cen1(k, 1), cen1(kk, 1)], [cen1(k, 2), cen1(kk, 2)]);
              
            %end;
    end;
end;

% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:n 
    ss=['station',num2str(k)];   
    TH=[];
    P1=cen1(k,:);
    neighbor=find( (dij(k,:) > 0) & (dij(k,:) <=700));  
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

% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xij=zeros(n);
for k=1:n 
    ss=['station',num2str(k)];   
    k1=eval([ss,'.neighborsNo']);   
    xij(k,k1)=1;
end;       

% Generate random traffic demands for the requests
m_lambda = 150; % Average traffic amount per request
requests = cell(n, 1);

for i = 1:n
    k = poissrnd(1); % Number of requests originating from node i
    destinations = setdiff(1:n, i); % Pick destinations from other nodes
    dest_indices = randsample(destinations, k, true);
    traffic_demands = normrnd(m_lambda, sqrt(0.5 * m_lambda), k, 1);
    %(squr(0.5*m_lambda) Standard deviation of the traffic demand between nodes
    requests{i} = [i * ones(k, 1), dest_indices', traffic_demands];
end

% Plot the nodes and requests
scatter(xc, yc,'filled');
hold on;
for i = 1:n
    for j = 1:size(requests{i}, 1)
        s = requests{i}(j, 1);
        d = requests{i}(j, 2);
        plot([xc(s), xc(d)], [yc(s), yc(d)], 'b');
    end
end
hold off;
title('Topology for Non-Splittable and Splittable Traffic');
xlabel('X');
ylabel('Y');
