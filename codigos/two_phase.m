% Escoamento oleo - agua (2-D) utilizando metodo IMPES
% Saturacao resolvida pelo metodo Upwind de primeira ordem
% Estos codigos ja estao em GitHub
clear all
clc
format short
global coord centelem elem esurn1 esurn2 nsurn1 nsurn2 bedge inedge ...
    normals esureface1 esureface2 esurefull1 esurefull2 elemarea dens ...
    visc satlimit pormap bcflag courant totaltime filepath foldername;
%%========================================================================%

[coord,centelem,elem,esurn1,esurn2,nsurn1,nsurn2,bedge,inedge,normals,...
    esureface1,esureface2,esurefull1,esurefull2,elemarea,dens,visc,...
    satlimit,pormap,bcflag,courant,totaltime,filepath,foldername,kmap,...
    wells] = preprocessor;
%[auxcoord]=distortedramd;
%% s� em malha com furo
% x=bedge(129:144,1);
% y=bedge(129:144,2);
% bedge(129:144,1)=y;
% bedge(129:144,2)=x;
%% s� em malha "Tipo1malha1"
% x=bedge(:,1);
% y=bedge(:,2);
% bedge(:,1)=y;
% bedge(:,2)=x;
% x1=elem(:,1);
% x2=elem(:,3);
% elem(:,1)=x2;
% elem(:,3)=x1;
%% no problema de Queiroz
% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
% x=bedge(71:78,1);
% y=bedge(71:78,2);
% bedge(71:78,1)=y;
% bedge(71:78,2)=x;
% bedge(71:78,4:5)=102; % 36x36
% bcflag(2,1)=102;
% bcflag(2,2)=2;
%-----------------------------
% unstructured mesh com furo reordenando o sentido da fronterira no
% contorno interior
%x=bedge(135:150,1);
%y=bedge(135:150,2);
%bedge(135:150,1)=y;
%bedge(135:150,2)=x;
%bedge(135:150,4:5)=102; % 36x36
%bcflag(2,1)=102;
%bcflag(2,2)=2;
%% Modifica��o Malha Kershaw
%bedge(:,4:5)=101;
%----------------------------
%% tratamento malha Hermeline
%bedge(:,4:5)=101;
% malha 16x16
% x=bedge(16:24,1);
% y=bedge(16:24,2);
% bedge(16:24,1)=y;
% bedge(16:24,2)=x;
% malha 32x32
% x=bedge(33:48,1);
% y=bedge(33:48,2);
% bedge(33:48,1)=y;
% bedge(33:48,2)=x;
% malha 64x64
% x=bedge(64:96,1);
% y=bedge(64:96,2);
% bedge(64:96,1)=y;
% bedge(64:96,2)=x;
% malha 128x128
%x=bedge(128:192,1);
%y=bedge(128:192,2);
%bedge(128:192,1)=y;
%bedge(128:192,2)=x;
%%
tic
%% escolha o m�todo de interpola��o 
% 1-->LPEW1 Gao e Wu 2010
% 2-->LPEW2 Gao e Wu 2010
interptype=2;
% num�rodo de Courant
CFL=courant;
%% exponente das permeabilidades relativas cuidado...! que varia com o benchmark
nw=2;
no=2;
%% Digite 
% monofasico ---> quando deseje rodar um problema de escoamento monof�sico ou
% bifasico   ---> quando deseja rodar um problema de ecoamento bif�sico 
simu='bifasico';
%% escolha o tipo de erro discreto que deseja usar
% erromethod1 ---> erro utilizado por Gao e Wu 2010
% erromethod2 --->  ''     ''     por Lipnikov et al 2010
% erromethod3 --->  ''     ''     por Eigestad et al 2005
% erromethod4 --->  ''     ''     por Shen e Yuan 2015
erromethod='erromethod1';
%% Defina o tipo de solver de press�o
% mpfad--> m�todo linear dos volumes finitos baseado no Gao e Wu 2010
% tpfa ---> m�todo Linear dos volumes finito 
% nlfvLPEW --> m�todo n�o linear dos volumes finitos baseado no artigo  (Sheng e Yuan, 2011)
% nlfvDMPSY --> m�todo n�o linear que preserva DMP baseado no artigo (Gao e Wu, 2013) e (Sheng e Yuan, 20...)
% lfvLPEW --> m�todo linear based no m�todo n�o linear usando LPEW;
% lfvHP --> m�todo linear baseado no m�todo n�o linear usando pontos
% harm�nicos 
pmetodo='mpfad';
%pmetodo='nlfvLPEW';
%% metodo de intera��o: iterpicard, iternewton, iterbroyden, itersecant,
%iterfreejacobian,iterdiscretnewton, JFNK
%iteration='iterdiscretnewton';
%iteration='iterpicard';
%iteration='iterbroyden';
iteration='iterpicard';
%iteration='iterhybrid';
%% defina o ponto de interpola��o
interpol='LPEW2';
%% Defina o tipo de solver de satura��o
%1. escreve ['FOU'] se deseja executar FOU method
%2. escreve 'HOFV-E' se deseja executar m�todo de alta por simple
%extrapola��o
%3. escreve 'HOMFV' se deseja executar m�todo de alta ordem com limitador
%de n�, na esque�a de definir upsilon (default � 0.2 como indica o artigo)
%2. escreve ['GOE'] se deseja executar o m�todo GOE-free;
smetodo='FOU';
order=1;
%% ordem de Runge Kutta
tordem=1;
% uso especial quando esta ativado o limitador de Woodfield, forne�a um
% valor no intervalo [0,1]
upsilon=0.2;
% Quando
% kapp=0--> m�todo de Fromm
% kappa=1--> m�todo de diferen�as centradas com tr�s pontos
% kappa=-1--> m�todo de pondera��o a montante de segunda ordem
% kappa=1/3--> m�todo de pondera��o a montante de terceira ordem
kappa=1/3;
% digite segundo o benchmark
% benchmark
% nikitin
% lamine
% durlofsky
% shuec
% buckley
benchmark='durlofsky';
% nome do arquivo � unico para cada exemplo
namefile='Report_Production_Mesh_lamine_LFVHP.dat';
% escreve sobre o arquivo criado
fid3 = fopen(namefile,'w');

%% adequa��o dos flags
%nflag= calflag(pmetodo);
%% este flag s� use quando o problema � Buckley-Leverett com fluxo imposto
% na face com flag 202, por�m a satura��o ser� imposto na face
auxflag=202;
%% adequa��o das permeabilidades segundo o bechmark
[elem,kmap,normKmap,solanal,bedge,fonte,velanal]=adequapermeab(benchmark, kmap,elem,bedge);
sat=ones(size(elem,1),1);
postprocessor(solanal,sat,1);

%% pre-precessador-satura��o
[N,F,V,weightLS,esuel1,esuel_coord,A,bound,S_old,S_cont]=presaturation(wells);

%% pre-m�todo-multidimensional
%[elemphant,elembedge,face_bedge,peso,nu]=premultidimensional(N);

%% pre-metodo-nao-linear
[pointarmonic,parameter,gamma,p_old,tol,nit,er,nflagface,nflagno,...
    weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface]=preNLFV(kmap,N,pmetodo,benchmark,interpol,bedge);
nflag=nflagno;
%[aroundface]=aroundfacelement(F,pointarmonic);
%% IMPES
% inicializando as vari�veis
cont = 1;
vpi_old = 0;
VPI(1) = 0;
cumulateoil(1) = 0;
oilrecovery(1) = 1;
watercut(1) = 0;
vpi = totaltime(2);
t_old = 0;
step=0;
time2=0;

while vpi_old < vpi
%while t_old<totaltime(2)
    %% calculo das mobilidades 
    [mobility] = mobilityface(S_old,nw,no,auxflag,S_cont,simu);
    
    %%  Calculo da Press�o Implicita pelos m�todo lineares e n�o-lineares
    [p,errorelativo,flowrate,flowresult]=solverpressure(kmap,nflagface,nflagno,fonte, tol,...
        nit,p_old,mobility,gamma,wells,S_old,V,nw,no,N,parameter,...
        pmetodo,auxflag,interptype,Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,benchmark,iteration,nflag);
        
    %% quando a velocidade � unitario
    %[p,influx,q] = getknownflowrate;
    %% calculo do fluxo fracional em cada elemento da malha
    f_elem = fractionalflow(S_old,nw,no);
    
    %% caculo do passo de tempo
    d_t=steptime(f_elem,S_old,flowrate,CFL,S_cont,auxflag,nw,no,order);
    
    %% calculo da satura��o explicito
    S_old=solversaturation(S_old,flowrate,d_t,esuel1,wells,flowresult,...
        f_elem,S_cont,weightLS,bound,smetodo,tordem,...
        upsilon,kappa,nw,no,auxflag);
        
    %% reporte de produ��o
    [VPI,oilrecovery,cumulateoil,watercut]=reportproduction(S_old,...
        wells,f_elem,cont,VPI,oilrecovery,cumulateoil,watercut,flowresult,d_t);
    
    %% calculo o passo de vpi ou passo de tempo
    %t_old=VPI
    
    vpi_old=VPI(cont)
   
    %visualiza��o
    if vpi_old >= time2
       step = 1000*time2;
       postprocessor(p,S_old,step)
      time2 = time2 + 0.01;
    end
    % armezena a produ��o em arquivo .dat
    fprintf(fid3,'%13.3e %13.3e %13.3e %13.3e \n',VPI(cont), ...
       oilrecovery(cont),cumulateoil(cont),watercut(cont));
    vpi_old=1e+40
   %postprocessor(p,S_old,cont+1)
   cont=cont+1
end
fclose(fid3);
%% calculo do erro
[erropressure,errovelocity]=errorateconv(solanal, p, velanal,flowrate,erromethod)
%% calculo das press�es m�ximas e minimas
panalmax= max(p)
panalmin= min(p)

toc

% if max(wells)==0
% plotBL(S_old,h)
% end
% S_old;

