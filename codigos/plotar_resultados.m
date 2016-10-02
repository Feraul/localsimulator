%% Resultados Multi-D Kozdon Explicito via Marcio
fid = fopen('ProdutionReport_reference_MPFA_O_FOU.dat');
b1 = fscanf(fid,'%g %g %g %g %g',[5 inf]); % It has two rows now.
b1=b1';
fclose(fid)

fid = fopen('Report_Production_Mesh_lamine_NLFV_FOU.dat');
b2 = fscanf(fid,'%g %g %g %g ',[4 inf]); % It has two rows now.
b2=b2';
fclose(fid)

fid = fopen('Report_Production_Mesh_lamine_NLFV_HOMFV.dat');
b3 = fscanf(fid,'%g %g %g %g ',[4 inf]); % It has two rows now.
b3=b3';
fclose(fid)

fid = fopen('Report_Production_Mesh_lamine_TPFA_FOU.dat');
b4 = fscanf(fid,'%g %g %g %g ',[4 inf]); % It has two rows now.
b4=b4';
fclose(fid)


%%
figure(1)

plot(b2(:,1),b2(:,2),'b','LineWidth',1.5)
hold on

plot(b3(:,1),b3(:,2),'r','LineWidth',1.5)
hold on

plot(b4(:,1),b4(:,2),'g','LineWidth',1.5)
hold on

plot(b1(:,1),b1(:,2),'k','LineWidth',1.5)

xlabel('Pore Volumes Injected')
ylabel('Oil Recovery')
legend('NLFV/FOU','NLFV/HOMFV','TPFA/FOU','Reference')
%legend('MPFA-D/FOU method','MPFA-D/MUSCL-LWf com -1','HO-Durlofsky','MPFA-D/MUSCL-LWf com 0.3333','Reference')

grid

%%
figure(2)

plot(b2(:,1),b2(:,3),'b','LineWidth',1.5)
hold on

plot(b3(:,1),b3(:,3),'r','LineWidth',1.5)
hold on

plot(b4(:,1),b4(:,3),'g','LineWidth',1.5)
hold on

plot(b1(:,1),b1(:,3),'k','LineWidth',1.5)

xlabel('Pore Volumes Injected')
ylabel('Cumulative Oil')
legend('NLFV/FOU','NLFV/HOMFV','TPFA/FOU','Reference')
%legend('MPFA-D/FOU method','MPFA-D/MUSCL-LWf com -1','HO-Durlofsky','MPFA-D/MUSCL-LWf com 0.333','Reference')

grid

% figure(3)
% plot(b2(:,1),b2(:,4),'b','LineWidth',1.5)
% hold on
%
% plot(b4(:,1),b4(:,4),'g','LineWidth',1.5)
% hold on
%
%  plot(b5(:,1),b5(:,4),'m','LineWidth',1.5)
%  hold on
%
%  plot(b3(:,1),b3(:,4),'r','LineWidth',1.5)
%  hold on
%
% plot(b1(:,1),b1(:,4),'k','LineWidth',1.5)
%
% xlabel('Pore Volumes Injected')
% ylabel('Water Cut')
% legend('MPFA-D/FOU method','MPFA-D/MUSCL-LWf','MPFA-D/MUSCL-NLWf','HO-Durlofsky','Reference')
% %legend('MPFA-D/FOU method','MPFA-D/MUSCL-¨LWf com -1','HO-Durlofsky','MPFA-D/MUSCL-¨LWf com 0.3333','Reference')
%
% grid
%
% figure(4)
% plot(b2(:,1),b2(:,5),'b','LineWidth',1.5)
% hold on
%
% plot(b4(:,1),b4(:,5),'g','LineWidth',1.5)
% hold on
%
% plot(b5(:,1),b5(:,5),'m','LineWidth',1.5)
% hold on
%
% plot(b3(:,1),b3(:,5),'r','LineWidth',1.5)
% hold on
%
% plot(b1(:,1),b1(:,5),'k','LineWidth',1.5)
%
% xlabel('Pore Volumes Injected')
% ylabel('Water Cut')
% legend('MPFA-D/FOU method','MPFA-D/MUSCL-NLWf','MPFA-D/MUSCL-LWf','HO-Durlofsky','Reference')
%
% grid
