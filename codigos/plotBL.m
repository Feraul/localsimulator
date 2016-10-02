function []=plotBL(S_old,h)
global bedge centelem 
s=h;
x=zeros(s,1);
k=1;
for l=1:size(bedge,1)
    
    if coord(bedge(l,1),2)==0 && bedge(l,5)==201
        x(k)=centelem(bedge(l,3),1);
        S(k)=S_old(bedge(l,3));
        k=k+1;
    end
    
end
c(1)=0;
c(2:s+1)=x(1:s);
Satur(1)=S_cont;
Satur(2:s+1)=S(1:s);
plot(c(1:s+1),Satur,'k-','LineWidth',1.5)
hold on
grid
BL=[0.00000 0.900000
    0.013750	0.882222
    0.027500	0.863175
    0.056250	0.831429
    0.086250	0.803492
    0.123750	0.774286
    0.167500	0.743810
    0.218750	0.713333
    0.300000	0.667619
    0.300000	0.100000
    1.000000	0.100000];

plot(BL(:,1),BL(:,2),'k-','LineWidth',1)
hold on
%legend('First order', 'Reference','High order')
xlabel('x-axes')
ylabel('Sw')
%fclose(fid3);
end