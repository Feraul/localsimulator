function [p,step,errorelativo,flowrate,flowresult]=iterpicard(M_old,RHS_old,nitpicard,tolpicard,kmap,...
    parameter,metodoP,auxflag,w,s,nflagface,fonte,p_old,gamma,nflagno,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded)

%% calculo do residuo Inicial
R0=norm(M_old*p_old-RHS_old);

%% inicializando dados para itera��o Picard
step=0;
er=1;
while (tolpicard<er || tolpicard==er) && (step<nitpicard)
    %% atualiza itera��es
    step=step+1
    %% calculo da press�o utilizando o precondicionador
    
    %[L,U] = ilu(M_old,struct('type','crout','droptol',1e-6));
    %M=L*U;
    
    %u=(M_old*inv(M))\RHS_old;
    
    %p_new=M\u;
    %% calculo da press�o pelo m�todo direto
    %p_new=inv(M_old)*RHS_old;  % invers�o com pivotamento
    
    p_new=M_old\RHS_old;  % invers�o sem pivotamento
    
    %% calculo da press�o usando m�todo iterativo + precondicionador ILU
    %[L,U]=ilu(M_old,struct('type','ilutp','droptol',1e-6));
    %M=L*U;
    %% GMRES do artigo SIAM
    %[p_new, error, iter, flag] = gmresSIAM( M_old, p_old, RHS_old, M, 10, 1000, 1e-10 );
    
    %% GMRES do MATLAB
    %[p_new,flag]=gmres(M_old,RHS_old,2,1e-6,100,L,U); % gmres Matlab
    %flag
    
    %% plotagem no visit
     S=ones(size(p_new,1),1);
     postprocessor(p_new,S,step)
    %% Interpola��o das press�es na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);
                                      
    %% Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflagface,nflagno...
        ,parameter,kmap,fonte,metodoP,w,s,benchmark,weightDMP,auxface,wells,mobility,Hesq, Kde, Kn, Kt, Ded);
    
    %% Calculo do residuo
    R = norm(M_new*p_new - RHS_new);
    
    %% calculo do erro
    A=logical(R0 ~= 0.0);
    B=logical(R0 == 0.0);
    er=A*abs(R/R0)+B*0
    
    errorelativo(step)=er;
    %% atualizar
    M_old=M_new;
    RHS_old=RHS_new;
    
end
p=M_old\RHS_old;
pinterp=pressureinterp(p,nflagface,nflagno,w,s,auxflag,metodoP,parameter,weightDMP);

if strcmp(metodoP,'nlfvDMPSY')
    % implementa��o do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(p, pinterp, parameter,nflagface,kmap,gamma,weightDMP,mobility);
else
    % implementa��o do fluxo NLFV
    [flowrate,flowresult]=flowrateNLFV(p, pinterp, parameter,mobility);
end
residuo=er;
niteracoes=step;

name = metodoP;
X = sprintf('Calculo o campo de press�o pelo m�todo: %s ',name);
disp(X)

x=['Erro:',num2str(residuo)];
disp(x);
y=['N�mero de itera��es:',num2str(niteracoes)];
disp(y);

end