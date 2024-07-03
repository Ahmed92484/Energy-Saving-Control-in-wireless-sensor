clear   
close all
warning off
clc
format short g   
Band=500;   
counter=0;   
load centers_n_15_2
string_table=[];   
tic
% Energy Constants
%%%%%%%%%%%%%%%%%%%
Eamp =100*10^(-12);   
packet_size=2000;   

neighbor_limit=150;  
Small_Stations(cen1,num1,limits,neighbor_limit);  

% Distance Matrix
dij=zeros(n,n);   

for k=1:n   
    for kk=1:n   
        dij(k,kk)=sqrt((cen1(k,1)-cen1(kk,1))^2+(cen1(k,2)-cen1(kk,2))^2);    
    end   
end   

% Cost Matrix
%%%%%%%%%%%%%%
Cij=round(dij/10);     


% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:n 
    ss=['station',num2str(k)];   
    TH=[];
    P1=cen1(k,:);
    neighbor=find( (dij(k,:) > 0) & (dij(k,:) <= neighbor_limit));  
    if (length(neighbor) >=1)
        for kk=1:length(neighbor)
            P2=cen1(neighbor(kk),:);   
            theta=theta1(P1,P2);
            TH(kk)=theta;  
        end   
        eval([ss,'.neighborsNo = neighbor;']);   
        eval([ss,'.neighborsAng = TH;']);  
    else
        eval([ss,'.neighborsNo = [];']);   
        eval([ss,'.neighborsAng = [];']);  
    end    
end  

% Calculate the stations within range to each one
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xij=zeros(n);
for k=1:n 
    ss=['station',num2str(k)];   
    k1=eval([ss,'.neighborsNo']);   
    xij(k,k1)=1;
end       

% assuming a value to alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=2;    
         
% Generate random numbers with mean lamdam and 
% variance lamdam/2
% Pick a request lamda_s_d
% 
%

lamdam=150;  
rr = lamdam + sqrt(lamdam/2).*randn(10000,1);

% var(r),mean(r)

% The traffic demand between node pair (s,d) lamda_sd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% random request generation (s is row, d is column)
sd=randi([0 1],n);   
for k=1:n   
    sd(k,k)=0;
end         

% Energy initialization
% each node with 0.5 Joul
%%%%%%%%%%%%%%%%%%%%%%%%
Energy_Vector=0.5*ones(1,n);  
for kk=1:n
eval(['station',num2str(kk),'.energy=0.5;'])
end   

% consider number of requests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nr=size(reqmat,1);   

% Take part of the requests
%%%%%%%%%%%%%%%%%%%%%%%%%%%
 nr=5;    
for counter1=1:nr   

lamda_sd=reqmat(counter1,1);   
s=reqmat(counter1,2);   
d=reqmat(counter1,3);   
V=1:n;    
OV=V;
i_s=find(V~=s);
V=V(i_s);
i_d=find(V~=d);
V=V(i_d);
rr=[];   
% c1=0;    
% Consider number of hops
%%%%%%%%%%%%%%%%%%%%%%%%%
delta=floor((n-1)/3);   
for i=1:delta   
    cc=combnk(V,i);
    for ii=1:size(cc,1)
        pp=perms(cc(ii,:));    
        Route_Path=[pp,zeros(size(pp,1),n-size(pp,2)-2)];   
        rr=[rr;Route_Path];
    end     
end
% Remove the zeros and add the source and destination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
route_number=0;   
% The first route from source to destination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rr=[[s d,zeros(1,size(rr,2)-2)];rr];     
route_distance=[];    
for i=1:size(rr,1)   
    ii=find(rr(i,:) ~= 0);
    if i>1
        route1=[s rr(i,ii) d];   
    else
        route1=rr(i,ii);   
    end      
    route_distance1=0; 
    Valid_route=1;   
    for vv=1:length(route1)-1
        route_distance1=route_distance1+dij(route1(vv),route1(vv+1));   
        vv1=route1(vv);   
        check_neighbor=eval(['find(station',num2str(vv1), ...
            '.neighborsNo == route1(',num2str(vv+1),'));']);           
        if length(check_neighbor) == 1
            check_neighbor1=1;  
        else
            check_neighbor1=0;  
        end
        Valid_route=Valid_route*check_neighbor1;   
    end
    % Check for being neighbors
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    if Valid_route
        route_number=route_number+1;    
        route_distance(route_number)=route_distance1;           
        eval(['route_number_',num2str(route_number), ...
            '_n_eq_',num2str(n),'_s_eq_',num2str(s), ...
            '_d_eq_',num2str(d),'.route=route1;'])
        eval(['route_number_',num2str(route_number), ...
            '_n_eq_',num2str(n),'_s_eq_',num2str(s), ...
            '_d_eq_',num2str(d),'.distance=route_distance1;']);                       
    end
end   

% Select the optimum route
%%%%%%%%%%%%%%%%%%%%%%%%%%%
selected_route_number=find(route_distance == min(route_distance));   
selected_route_number=selected_route_number(1);      
SRN=selected_route_number;   
eval(['opt_route_n_eq_',num2str(n),'_s_eq_', ...
    num2str(s),'_d_eq_',num2str(d), ...
    '=route_number_',num2str(SRN), ...
    '_n_eq_',num2str(n),'_s_eq_',num2str(s), ...
    '_d_eq_',num2str(d),'.route;'])
 
eval(['Route_Path=opt_route_n_eq_',num2str(n),'_s_eq_', ...
    num2str(s),'_d_eq_',num2str(d)]);

% Energy update for each node in the route
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First node in the route
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Et = Eamp*d^alpha* b (2)
% Et in Equation (2) is the energy sent by a node to transfer
% a series of bit. Eamp = 100 pj/bit/m2, d is the distance for

Et=Eamp*dij(Route_Path(1),Route_Path(2))^2*packet_size;   
eval(['station',num2str(Route_Path(1)),'.energy=', ...
    'station',num2str(Route_Path(1)),'.energy-Et;']);

% Last node in the route
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Er = Eelec*b (3)
% In Equation (3), the energy spent by a node to receive a
% packet is given Eelec = 50nj/bit
Eelec = 50*10^-9;  
Er = Eelec*packet_size;   
eval(['station',num2str(Route_Path(end)),'.energy=', ...
    'station',num2str(Route_Path(end)),'.energy-Et;']);

if length(Route_Path)>2
    for count2=2:length(Route_Path)-1
        Et=Eamp*dij(Route_Path(count2),Route_Path(count2+1))^2*packet_size;   
        Er = Eelec*packet_size;  
        Etotal=Et+Er;   
        eval(['station',num2str(Route_Path(count2)),'.energy=', ...
            'station',num2str(Route_Path(count2)),'.energy-Etotal;']);
    end   
end         
TotalCost=0;  
for count3=1:length(Route_Path)-1   
    TotalCost=TotalCost+Cij(Route_Path(count3),Route_Path(count3+1));   
end      

Total_Cost1(counter1)=TotalCost;   

Energy_Vector1=[];    
for cc1=1:n
eval(['Energy_Vector1=[Energy_Vector1 station',num2str(cc1),'.energy];']);
end    

Total_Consumed_Accumulated1(counter1)=sum(Energy_Vector-Energy_Vector1);        
 
path_string=[];  
for k=1:length(Route_Path)-1
    path_string=[path_string,num2str(Route_Path(k)),'-->' ];
end
path_string=[path_string,num2str(Route_Path(end))];    
counter=counter+1;   
table1(counter,1)=n;   
table1(counter,2)=s;   
table1(counter,3)=d;   
table1(counter,4)=lamda_sd;   
table1(counter,5)=Band;   
string_table=strvcat(string_table,path_string);           
% figstr=['fig',num2str(counter),'_s_',num2str(s),'_d_',num2str(d)];  
% saveas(gcf,figstr,'jpg')    
Tend1(counter1)=toc;   
end        

figure;
h=plot(1:nr,Total_Consumed_Accumulated1,'LineWidth',2);
xlabel('Request Number');
ylabel('Consumed Energy in Joul');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write the data to a text file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen('route_selection4.txt','w');   
fprintf(fid,'%70s\n\n\n',datestr(now,0));   
fprintf(fid,'\t\t Ad Hoc, Power Saving Analysis \n');  
fprintf(fid,'=================== Data Analysis ====================\n');  
fprintf(fid,'==========================================================\n\n'); 
fprintf(fid,'n\t s\t d\t lamda_sd\t Bandwidth \t\t  route\r\n');   
fprintf(fid,'--- --- --- ----------\t -------- -----------------------\t \r\n');      
for k=1:size(table1,1)  
   fprintf(fid,'%2d  %2d  %2d  %10.4f  %7d \t  %20s \r\n', ...
       table1(k,1),table1(k,2),table1(k,3),table1(k,4),table1(k,5), ...
            string_table(k,:));  
end        

fclose(fid);     

% save route_selection4_Data


