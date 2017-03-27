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
%% so em malha com furo
% x=bedge(129:144,1);
% y=bedge(129:144,2);
% bedge(129:144,1)=y;
% bedge(129:144,2)=x;
%% so em malha "Tipo1malha1"
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
%% Kershaw Mesh Modification
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
%% escolha o metodo de interpolacao
% 1-->LPEW1 Gao e Wu 2010
% 2-->LPEW2 Gao e Wu 2010
interptype=2;
% numero de Courant
CFL=courant;
%% exponente das permeabilidades relativas cuidado...! que varia com o benchmark
nw=2;
no=2;
%% Digite
% monofasico ---> quando deseje rodar um problema de escoamento monofï¿½sico ou
% bifasico   ---> quando deseja rodar um problema de ecoamento bifï¿½sico
simu='monofasico';
%% Tipo de Problema
% stationary   --> problema simple de diffusion sem tempo
% nonstationary --> problema de diffusion com passo de tempo
problemtype='nonstationary';
%% escolha o tipo de erro discreto que deseja usar
% erromethod1 ---> erro utilizado por Gao e Wu 2010
% erromethod2 --->  ''     ''     por Lipnikov et al 2010
% erromethod3 --->  ''     ''     por Eigestad et al 2005
% erromethod4 --->  ''     ''     por Shen e Yuan 2015
erromethod='erromethod1';
%% Defina o tipo de solver de pressao
% mpfad--> metodo linear dos volumes finitos baseado no Gao e Wu 2010
% tpfa ---> metodo Linear dos volumes finito
% nlfvLPEW --> metodo nao linear dos volumes finitos baseado no artigo  (Sheng e Yuan, 2011)
% nlfvDMPSY --> metodo nao linear que preserva DMP baseado no artigo (Gao e Wu, 2013) e (Sheng e Yuan, 20...)
% lfvLPEW --> metodo linear based no metodo nao linear usando LPEW;
% lfvHP --> metodo linear baseado no metodo nao linear usando pontos
% harmonicos
pmetodo='mpfad';
%pmetodo='nlfvLPEW';
%% metodo de interacao: iterpicard, iternewton, iterbroyden, itersecant,
%iterfreejacobian,iterdiscretnewton, JFNK
%iteration='iterdiscretnewton';
%iteration='iterpicard';
%iteration='iterbroyden';
iteration='iterpicard';
%iteration='iterhybrid';
%% defina o ponto de interpolacao
interpol='LPEW2';
%% Defina o tipo de solver de saturacao
%1. escreve ['FOU'] se deseja executar FOU method
%2. escreve 'HOFV-E' se deseja executar metodo de alta por simple
%extrapolacao
%3. escreve 'HOMFV' se deseja executar metodo de alta ordem com limitador
%de no, na esqueca de definir upsilon (default e 0.2 como indica o artigo)
%2. escreve ['GOE'] se deseja executar o metodo GOE-free;
smetodo='FOU';
order=1;
%% ordem de Runge Kutta
tordem=1;
% uso especial quando esta ativado o limitador de Woodfield, forneca um
% valor no intervalo [0,1]
upsilon=0.2;
% Quando
% kapp=0--> metodo de Fromm
% kappa=1--> metodo de diferencas centradas com tres pontos
% kappa=-1--> metodo de ponderacao a montante de segunda ordem
% kappa=1/3--> metodo de ponderacao a montante de terceira ordem
kappa=1/3;
% digite segundo o benchmark
% benchmark
% nikitin
% lamine
% durlofsky
% shuec
% buckley
benchmark='yuansheng2008';
% nome do arquivo  unico para cada exemplo
namefile='Report_Production_Mesh_lamine_LFVHP.dat';
% escreve sobre o arquivo criado
fid3 = fopen(namefile,'w');

%% adequacao dos flags
%nflag= calflag(pmetodo);
%% este flag se use quando o problema e Buckley-Leverett com fluxo imposto
% na face com flag 202, porm a saturacao sera imposto na face
auxflag=202;
%% adequacao das permeabilidades segundo o bechmark
[elem,kmap,normKmap,solanal,bedge,fonte,velanal]=adequapermeab(benchmark, kmap,elem,bedge);
sat=ones(size(elem,1),1);
%postprocessor(solanal,sat,1);

%% pre-precessador-saturacao
[N,F,V,weightLS,esuel1,esuel_coord,A,bound,S_old,S_cont]=presaturation(wells);

%% pre-metodo-multidimensional
%[elemphant,elembedge,face_bedge,peso,nu]=premultidimensional(N);

%% pre-metodo-nao-linear
% modificação para o problema "yuansheng2008"
elem(:,5)=1;
[pointarmonic,parameter,gamma,p_old,tol,nit,er,nflagface,nflagno,...
    weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface]=preNLFV(kmap,N,pmetodo,benchmark,bedge);
nflag=nflagno;
for ielem=1:size(centelem,1)
    
p_old(ielem,1)=sin(centelem(ielem,1)*pi)*sin(centelem(ielem,2)*pi);
end
%[aroundface]=aroundfacelement(F,pointarmonic);
%% IMPES
% inicializando as variaveis
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
p_anal=zeros(size(elem,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f_elem=0;

while t_old<totaltime(2)
    p_anal=exp(-2*(pi^2)*t_old)*sin(pi*centelem(:,1)).*sin(pi*centelem(:,2));
    %% calculo das mobilidades
    [mobility] = mobilityface(S_old,nw,no,auxflag,S_cont,simu);
    %%  Calculo da Pressao Implicita pelos matodo lineares e nao-lineares
    [~,~,flowrate,flowresult]=solverpressure(kmap,nflagface,nflagno,fonte, tol,...
        nit,p_old,mobility,gamma,wells,S_old,V,nw,no,N,parameter,...
        pmetodo,auxflag,interptype,Hesq, Kde, Kn, Kt, Ded,weightDMP,auxface,benchmark,iteration,nflag,problemtype);
    
    %% caculo do passo de tempo
    %d_t=steptime(f_elem,p_old,flowrate,CFL,S_cont,auxflag,nw,no,order,problemtype);
    d_t=1/size(elem,1);
    %% calculo da pressão explicito
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    RHS=zeros(size(elem,1),1);
    
    %% aproximação upwind
    for iface = 1 : size(inedge,1)
        
        lef = inedge(iface,3); % elemento a esquerda
        rel = inedge(iface,4); % elemento a direita
        ve_mais = (flowrate(iface+size(bedge,1)) + abs(flowrate(iface+size(bedge,1))))/2;
        ve_menos = (flowrate(iface+size(bedge,1)) - abs(flowrate(iface+size(bedge,1))))/2;
        
        RHS(rel) = RHS(rel) + ve_mais + ve_menos;
        RHS(lef)  = RHS(lef)  - ve_mais - ve_menos;
        
    end
    %% calculo dos RHS nos poços produtores
    % calculo dos contribuições nos poços
    
    for ifacont = 1:size(bedge,1)
        lef=bedge(ifacont,3);
        if bedge(ifacont,5)==auxflag || bedge(ifacont,5)==102
            
            RHS(lef) = RHS(lef) -  flowrate(ifacont);
        elseif bedge(ifacont,5)<200
            
            RHS(bedge(ifacont,3)) = RHS(bedge(ifacont,3)) - flowrate(ifacont);
        end
    end
    
    %% calculo da saturação
    for i = 1:size(elem,1)
        
        p_new(i,1) = p_old(i,1) - (d_t*RHS(i))/(elemarea(i));
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    t_old=t_old+d_t;
    %visualizacao
    %if vpi_old >= time2
    %   step = 1000*time2;
    %   postprocessor(p,S_old,step)
    %  time2 = time2 + 0.01;
    %end
    
    postprocessor(p_anal,p_new,cont+1)
    p_old=p_new;
    cont=cont+1
end
fclose(fid3);
%% calculo do erro
[erropressure,errovelocity]=errorateconv(solanal, p, velanal,flowrate,erromethod)
%% calculo das pressoes maximas e minimas
panalmax= max(p)
panalmin= min(p)

toc

% if max(wells)==0
% plotBL(S_old,h)
% end
% S_old;
