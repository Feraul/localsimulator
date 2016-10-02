function [p,step,errorelativo,flowrate,flowresult]=iterbroyden(M_old,RHS_old,p_old,tolnewton,kmap,parameter,...
    metodoP,auxflag,w,s,nflag,fonte,gamma,nflagno,benchmark,R0_old,p_old1,weightDMP,auxface)

%% inicializando dados para iteração Picard
R0=M_old*p_old-RHS_old;
er=1;
step=0;
neta=1;
%Br=eye(size(p_old,1),size(p_old,1));
% primeira aproximação vai pela matriz Jacobiana
[Br]=aproxmjacobian(R0,p_old1,p_old,nflag,w,s,metodoP,parameter,kmap,nflagno,benchmark,fonte,auxflag,gamma,weightDMP,auxface);

while tolnewton<er
    step=step+1;
    % calculo da pressão
    sk=-inv(Br)*(R0_old);
    p_new=p_old1+sk;
    % plotagem no visit
    postprocessor(p_new,step)
    % Interpolação das pressões na arestas (faces)
    [pinterp_new]=pressureinterp(p_new,nflag,w,s,auxflag,metodoP,parameter,weightDMP);
    
    % Calculo da matriz global
    [M_new,RHS_new]=globalmatrix(p_new,pinterp_new,gamma,nflag,nflagno...
        ,parameter,kmap,fonte,metodoP,w,benchmark,weightDMP,auxface);
    % calculo do residuo
    R0_new= M_new*p_new - RHS_new;
    yk=R0_new-R0_old;
    oldBr=Br;
    %----------------------------------
    % good Broyden update
    if sk'*inv(oldBr)*yk~=0
        % original the Broyden formule
        % "neta" eh introduzido para que o Br sempre seja matriz no-singular
        Br=oldBr+neta*((yk-oldBr*sk)*sk')/(sk'*sk);
        % the BFGS formule
        %Br=oldBr-((oldBr*sk*sk'*oldBr)/(sk'*oldBr*sk))+ (yk*yk')/(yk'*sk);
        % the DFP formule
        %Br=(I-(yk*sk')/(yk'*sk))*oldBr*(I-(yk*sk')/(yk'*sk))'+(yk*yk')/(yk'*sk);
        
    else
        Br=oldBr;
    end
    %     tau=sk'*inv(Br)*yk/(sk'*sk);
    %     if abs(tau)>0.1|| abs(tau)==0.1
    %         neta=1;
    %     elseif abs(tau)<0.1
    %         neta=(1-0.1*sign(tau))/(1-tau);
    %     end
    % calculo do erro
    A=logical(norm(R0) ~= 0.0);
    B=logical(norm(R0) == 0.0);
    er=A*abs(norm(R0_new)/norm(R0))+B*0
    
    errorelativo(step)=er;
    R0_old=R0_new;
    p_old1=p_new;
end
p=p_old1;
pinterp=pressureinterp(p,nflag,w,s,auxflag,metodoP,parameter,weightDMP);
if strcmp(metodoP,'nlfvDMPSY')
% implementação do fluxo NLFV-DMP
    [flowrate,flowresult]=flowrateNLFVDMP(p, pinterp, parameter,nflag,kmap,gamma,weightDMP);
else
% implementação do fluxo NLFV
    [flowrate,flowresult]=flowrateNLFV(p, pinterp, parameter);
end
end