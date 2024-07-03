function []=stations_random_node(cen1,num1)   
% function []=stations_random_node(cen1,num1)   
% cen1 is the center of the stations
% num1 is the numbers
figure('position',[100 100 850 600]);   
xxc=[];  
yyc=[];  
neighbor_xc=[];  
neighbor_yc=[];  
for k=1:length(num1)
    [xc yc]=circle([cen1(k,1) cen1(k,2)],120);    
    [xn yn]=circle([cen1(k,1) cen1(k,2)],900);        
    xxc=[xxc,xc'];   
    yyc=[yyc,yc'];   
    neighbor_xc=[neighbor_xc,xn'];  
    neighbor_yc=[neighbor_yc,yn'];      
end;    
region_size = 3500;
xxc=[xxc,neighbor_xc];   
yyc=[yyc,neighbor_yc]; 
plot(xxc,yyc),axis image,
axis([0 region_size 0 region_size]);  



hold on
for k=1:length(num1)
    gt=text(cen1(k,1)-100,cen1(k,2)-50,num2str(num1(k)));   
    set(gt,'FontSize',12,'FontWeight','bold')
end;


title('Stations Random Node ');
xlabel('Units in meters');  
ylabel('Units in meters');  