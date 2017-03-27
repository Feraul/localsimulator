function [pointarmonic,parameter,gamma,p_old,tol,nit,er,nflagface,...
    nflagno,weightDMP,Hesq,Kde,Kn,Kt,Ded,auxface]=preNLFV(kmap,N,metodoP,benchmark,bedge)
global elem

nflagface=0;
pointarmonic=0;
parameter=0;
auxface=0;
weightDMP=0;
Hesq=0;
Kde=0;
Kn=0;
Kt=0;
Ded=0;

if strcmp(metodoP,'nlfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPS(kmap); % urgente revisar
    %temos usado para muitos estes o seguinte rutina
    [parameter]=coefficientLPSangle(kmap); 
    % adequação dos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    %% calculo dos pesos
    %[w,s] = weightnlfv(kmap,N,interpol);
elseif strcmp(metodoP,'nlfvDMPSY')|| strcmp(metodoP,'lfvHP') || strcmp(metodoP,'nlfvDMPV1')
    %% faces alrededor de um elemento
    [facelement]=element_face;
    %% calculoa dos pontos armonicos
    [pointarmonic]=harmonicopoint(kmap,N,benchmark);
    %% calculo dos parametros ou constantes (ksi)
    % temos usado este parametro durante muito tempo em muitos testes
    [parameter,auxface]=coefficientPPSharmonicpoint(facelement,pointarmonic,kmap);
    % esta rutina estamos usando de 7/2/2016
    %[parameter]=coefficientPPSusingHP(kmap,facelement,pointarmonic); %para lfvHP
    % adequação dos flag de face de contorno
    nflagface= contflagface(benchmark,bedge);
    % adequação dos nos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    %% calculo dos pesos DMP
    [weightDMP]=weightnlfvDMP(kmap);
elseif strcmp(metodoP,'lfvLPEW')
    %% calculo dos parametros ou constantes (ksi)
    
    %[parameter]=coefficientLPS(kmap);
    [parameter]=coefficientLPSangle(kmap);
    % adequação dos flags de contorno
    nflagno= contflagno(benchmark,bedge);
    % calculo dos pesos DMP
    [weightDMP]=weightnlfvDMP(kmap);

elseif strcmp(metodoP,'mpfad')
    
    %% calculo das constantes fisicos-geometrico
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
     % adequação dos flags de contorno
    nflagno= contflagno(benchmark,bedge);
else
    %% calculo das constantes fisicos-geometrico
    [Hesq, Kde, Kn, Kt, Ded]=preMPFAD(kmap);
     % adequação dos flags de contorno
    nflagno= contflagno(benchmark);
    
end
%% dados inicialização métodos dos volumes finitos não linear
gamma=0.0; % este parametro esta no intervalo [0,1]
% inicializando a pressão
p_old=1e1*ones(size(elem,1),1);
tol=1e-12; % tolerancia
nit=20000; % iterações
% Picard Iterations
er=1;
end