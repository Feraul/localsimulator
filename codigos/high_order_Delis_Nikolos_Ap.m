function [S,limiter_aqui,limiter_outro] = high_order_Delis_Nikolos_Ap(S,VN,dt,volume,satlimit,inedge,...
    c,visc,elem,esuel,bedge,coord,pormap,wells,q,f_elem,vn,bcflag,bflux,S_cont,Sat_max_min,w,bound,FF,xM)

% c: baricentro da célula
% gradient in element
%g = grad_S(elem,c,S,esuel1,esuel2,local_volume,viz_b,bedge,bcflag,inedge);
[g]=grad(S,w,esuel,elem,bcflag,bedge,bound);
RHS = zeros(size(S,1),1);
limiter_outro=zeros(size(inedge,1),1);
limiter_aqui=zeros(size(inedge,1),1);
for i = 1:size(inedge,1)
    
    v1 = inedge(i,1);
    v2 = inedge(i,2);
    
    outro=inedge(i,4);    %indice do elemento a direita da face i
    aqui=inedge(i,3);    %indice do elemento a esquerda da face i
    por=pormap(1);
    
    %% Definie o K
    K=0.3;
    %%
 list_outro = esuel(outro,:);

 list_aqui=esuel(aqui,:);

    t=1;
    for y=1:size(list_outro,2)
        if list_outro(y)~=outro
            p(t)=S(list_outro(y));
            t=t+1;
        end
       
    end
    
    delt_max_outro=max(S(outro),max(p))-S(outro);
    delt_min_outro=min(S(outro),min(p))-S(outro);
    e_outro= K*volume(outro);
    for m=1:size(list_outro,2)
        rPM_outro=FF(m,:,outro)-c(outro,:);
        delt_outro=sign(dot(g(list_outro(m),:),rPM_outro))*(abs(dot(g(list_outro(m),:),rPM_outro))+ eps);
        
        if delt_outro>0
            alfa_outro(m)=(1/delt_outro)*(((delt_max_outro^2+e_outro^2)*delt_outro +2*delt_outro^2*delt_max_outro)/(delt_max_outro^2+2*delt_outro^2+ delt_outro*delt_max_outro+e_outro^2));
        elseif delt_outro<0
            alfa_outro(m)=(1/delt_outro)*(((delt_min_outro^2+e_outro^2)*delt_outro +2*delt_outro^2*delt_min_outro)/(delt_min_outro^2+2*delt_outro^2+ delt_outro*delt_min_outro+e_outro^2));
        elseif delt_outro==0
            alfa_outro(m)=1;
        end
    end
    
    limiter_outro(i)=min(alfa_outro);
    w=1;
    for h=1:size(list_aqui,2)
        if list_aqui(h)~=aqui
            s(w)=S(list_aqui(h));
            w=w+1;
        end
       
    end
    
    delt_max_aqui=max(S(aqui),max(s))-S(aqui);
    delt_min_aqui=min(S(aqui),min(s))-S(aqui);
    
    e_aqui= K*volume(aqui);
    
    for m=1:size(list_aqui,2)
        rPM_aqui= FF(m,:,aqui)-c(aqui,:);
        delt_aqui=sign(dot(g(list_aqui(m),:),rPM_aqui))*(abs(dot(g(list_aqui(m),:),rPM_aqui))+ eps);
        
        if delt_aqui>0
            alfa_aqui(m)=(1/delt_aqui)*(((delt_max_aqui^2+e_aqui^2)*delt_aqui +2*delt_aqui^2*delt_max_aqui)/(delt_max_aqui^2+2*delt_aqui^2+ delt_aqui*delt_max_aqui+e_aqui^2));
        elseif delt_aqui<0
            alfa_aqui(m)=(1/delt_aqui)*(((delt_min_aqui^2+e_aqui^2)*delt_aqui +2*delt_aqui^2*delt_min_aqui)/(delt_min_aqui^2+2*delt_aqui^2+ delt_aqui*delt_min_aqui+e_aqui^2));
        elseif delt_aqui==0
            alfa_aqui(m)=1;
        end
    end
    limiter_aqui(i)=min(alfa_aqui);
    
%     SM_L = S(aqui) + limiter_aqui(i)*dot(g(aqui,:), xM(i,:)-c(aqui,:));
%     SM_R = S(outro) +limiter_outro(i)*dot(g(outro,:), xM(i,:)-c(outro,:));
    
    %%
    if Sat_max_min(outro,1)-Sat_max_min(outro,2)<=1e-20
        phi_outro= 1;
    else
        gamma_outro= (S(outro)-Sat_max_min(outro,2))/(Sat_max_min(outro,1)-Sat_max_min(outro,2));
        if gamma_outro >=1 | gamma_outro<=0
            phi_outro=0;
        elseif 0.2<= gamma_outro && gamma_outro<=1-0.2
            phi_outro=1;
        elseif 0< gamma_outro && gamma_outro< 0.2
            phi_outro=gamma_outro/0.2;
        elseif 1-0.2 < gamma_outro && gamma_outro<1
            phi_outro=(1-gamma_outro)/0.2;
        end
    end
    
    if Sat_max_min(aqui,1)-Sat_max_min(aqui,2)<=1e-20
        phi_aqui= 1;
    else
        gamma_aqui= (S(aqui)-Sat_max_min(aqui,2))/(Sat_max_min(aqui,1)-Sat_max_min(aqui,2));
        if gamma_aqui >=1 | gamma_aqui<=0
            phi_aqui=0;
        elseif 0.2<= gamma_aqui && gamma_aqui<=1-0.2
            phi_aqui=1;
        elseif 0< gamma_aqui && gamma_aqui< 0.2
            phi_aqui=gamma_aqui/0.2;
        elseif 1-0.2 < gamma_aqui && gamma_aqui<1
            phi_aqui=(1-gamma_aqui)/0.2;
        end
    end
    
            SM_L = S(aqui) + phi_aqui*limiter_aqui(i)*dot(g(aqui,:), xM(i,:)-c(aqui,:));
            SM_R = S(outro) + phi_outro*limiter_outro(i)*dot(g(outro,:), xM(i,:)-c(outro,:));
    %%
    
    
    [fe_L] = flow_fract(visc,SM_L,satlimit); % f(SM_L)
    
    [fe_R] = flow_fract(visc,SM_R,satlimit); % f(SM_R)
    
    
    ve_mais = (VN(i) + abs(VN(i)))/2;
    ve_menos = (VN(i) - abs(VN(i)))/2;
    RHS(outro) = RHS(outro) + ve_mais*fe_L + ve_menos*fe_R;
    RHS(aqui)  = RHS(aqui)  - ve_mais*fe_L - ve_menos*fe_R;
end
for iw=1:size(elem,1)
      iwell=wells(iw,2);
    if iwell==2
        
        sink= f_elem(wells(iw,1))*q(wells(iw,1));
        RHS(wells(iw,1)) = RHS(wells(iw,1)) + sink;
    end
end
for i = 1:size(S,1)
    if S(i)~=1
        S(i,1) = S(i,1) + (dt*RHS(i))/(por*volume(i));
    end
end

end