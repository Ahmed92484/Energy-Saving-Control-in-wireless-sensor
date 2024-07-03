clear all
close all
clc
format compact
%%%%%%%%%%%%%%%%%%%
% number of nodes
%%%%%%%%%%%%%%%%%%%
% n=5;   
n=input('Number of Nodes \n');
%%%%%%%%%%%%%%%%%%%%%
V=1:n;   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assume a value for the source
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %s=1;
s=input('Source \n');
d=input('Destination \n');
%station_distanse=zeros()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assume a value for the destination
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% d=4;
xxx=[];
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
for i=1:size(rr,1)
    ii=find(rr(i,:) ~= 0);
    route_number=route_number+1;
    route=[s rr(i,ii) d];   
    eval(['route_number_',num2str(route_number),'=route'])
end;   
eval(['number_of_routes_for_n_eq_',num2str(n),'=route_number'])


