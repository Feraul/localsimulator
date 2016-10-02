function [S] = high_order_Delis_Nikolos(S,VN,dt,volume,satlimit,inedge,...
    c,visc,elem,xD,QI,QJ,bedge,coord,pormap,wells,q,f_elem,vn,bcflag,bflux,S_cont,Sat_max_min,w,bound,FF,esuel,xM)

% c: baricentro da célula
% gradient in element
%g = grad_S(elem,c,S,esuel1,esuel2,local_volume,viz_b,bedge,bcflag,inedge);
[g]=grad(S,w,esuel,elem,bcflag,bedge,bound);

% vetor termo de fonte
RHS = zeros(size(S,1),1);
for i = 1:size(inedge,1)
      
    v1 = inedge(i,1);
    v2 = inedge(i,2);
    
    rel=inedge(i,4);    %indice do elemento a direita da face i
    lef=inedge(i,3);    %indice do elemento a esquerda da face i
    por=pormap(1);
    
    S_aqui = S(lef);
    S_outro = S(rel);
    
    aqui=lef;
    outro=rel;
    
    
    rDM = (xM(i,:)-xD(i,:))';
    %%  calculo do gradiente Upwind
    
    
    gradS_cent_aqui = ( S_outro - S_aqui);
    
    r_aqui  =(c(outro,:)-c(aqui,:))';
    
    r_aqui_D = xD(i,:)-c(aqui,:);
    
    norma_aqui=(norm(r_aqui_D)/norm(r_aqui));
    
    prod_int2=dot(g(aqui,:)',r_aqui);
    
    gradS_upwind_aqui  = 2*prod_int2 - gradS_cent_aqui; % OK
    
    
    %   Saturation on point D
    if gradS_upwind_aqui*gradS_cent_aqui>0
        
        r=gradS_upwind_aqui/gradS_cent_aqui;
        
        limit_aqui=min((r^2+r)/(r^2+1), 1/(0.95*r));
        
        %limit_aqui(i)=min(r,1);
    else
        limit_aqui=0;
    end
    
    SD_aqui  =  S_aqui +  norma_aqui*limit_aqui*gradS_cent_aqui;
    
        
    %%  calculo do gradiente Upwind
    
    
    gradS_cent_outro = ( S_outro - S_aqui);
    
    r_outro =(c(aqui,:)-c(outro,:))';
    
    r_D_outro = c(outro,:)-xD(i,:);
    
    norma_outro=(norm(r_D_outro)/norm(r_outro));
    
    prod_int1=dot(g(outro,:)',r_outro);
    
    gradS_upwind_outro = 2*prod_int1 - gradS_cent_outro; % OK
    
    
    if gradS_upwind_outro*gradS_cent_outro>0
        
        rr=gradS_upwind_outro/gradS_cent_outro;
        
        limit_outro=min((rr^2+rr)/(rr^2+1),1/(rr*0.95));
        
        %limit_outro(i)=min(rr,1);
    else
        limit_outro=0;
    end
    SD_outro =  S_outro - norma_outro*limit_outro*gradS_cent_outro;
    
    
    %%  Saturation on point D
    
    %     limit_aqui(i)=limiter(gradS_upwind_aqui, gradS_cent_aqui); % OK
    %     limit_outro(i)=limiter(gradS_upwind_outro, gradS_cent_outro); % OK
    %
    %     SD_aqui  =  S_aqui +  norma_aqui*limit_aqui(i);
    %
    %     SD_outro =  S_outro + norma_outro*limit_outro(i);
    
    %% saturação no ponto M (medio na face)
    
    
    % Saturação sobre o ponto M
    
    % "M" é o ponto medio
    
    % "Q" é o ponto tal que o angulo entre rDM e rIQ é minimo
    
    % pk2 = pl2+l2k2 ---> l2k2 = pk2 - pl2 (veja figura 3)
    %if norm(rDM)<1e-16
        SM_aqui= SD_aqui ;
        
        
        SM_outro= SD_outro;
%     else
%         
%         rQI_proj = dot((c(QI(i),:)-c(aqui,:))',rDM)/(norm(rDM))^2*rDM -(c(QI(i),:)-c(aqui,:))';
%         
%         
%         S_QI = S(QI(i)) + dot(rQI_proj,g(QI(i),:)');
%         
%         % qm2 = r2m2 + pr2 ---> r2m2 = qm2 - pr2 (veja figura 3)
%         
%         rQJ_proj = dot((c(QJ(i),:)-c(outro,:))',rDM)/(norm(rDM))^2*rDM -...
%             (c(QJ(i),:)-c(outro,:))';
%         
%         
%         S_QJ = S(QJ(i)) + dot(rQJ_proj,g(QJ(i),:)');
%         
%         gradS_cent_Jele_rJQJ = S_QJ - S_outro;
%         
%         gradS_cent_Iele_rIQI = S_QI - S_aqui;
%         
%         rIQI = dot((c(QI(i),:)-c(aqui,:))',rDM)/(norm(rDM))^2*rDM;
%         rJQJ = dot((c(QJ(i),:)-c(outro,:))',rDM)/(norm(rDM))^2*rDM;
% %         
% %         % verifica se o sinal dos gradientes são opostos a sinal da velocidade
% %         
%         gradS_upwind_Jele_rJQJ = 2*dot(g(outro,:)',rJQJ) - gradS_cent_Jele_rJQJ;
%         gradS_upwind_Iele_rIQI = 2*dot(g(aqui,:)',rIQI) - gradS_cent_Iele_rIQI;
% %         %%
% %         % Correção da Saturação sobre o ponto M
% %         %         limit_aqui1(aqui)= limiter(gradS_upwind_Iele_rIQI,gradS_cent_Iele_rIQI);
% %         %         limit_outro1(outro)=limiter(gradS_upwind_Jele_rJQJ,gradS_cent_Jele_rJQJ);
% %         %
% %         %         SM_aqui= SD_aqui + (norm(rDM)/norm(rIQI))*limit_aqui1(i);
% %         %
% %         %
% %         %         SM_outro= SD_outro + (norm(rDM)/norm(rJQJ))*limit_outro1(i);
% %         
% %         %%
%         if gradS_cent_Iele_rIQI*gradS_upwind_Iele_rIQI>0
%             
%             r1=gradS_upwind_Iele_rIQI/gradS_cent_Iele_rIQI;
%             
%             limit_aqui1=(r1^2+r1)/(r1^2+1);
%             
%             %limit_aqui1(i)=max(0,min(r1(i),1));
%         else
%             
%             limit_aqui1=0;
%         end
%         
%         if gradS_upwind_Jele_rJQJ*gradS_cent_Jele_rJQJ>0
%             
%             rr1=gradS_upwind_Jele_rJQJ/gradS_cent_Jele_rJQJ;
%             
%             limit_outro1=(rr1^2+rr1)/(rr1^2+1);
%             
%             %limit_outro1(i)=max(0,min(rr1(i),1));
%         else
%             
%             limit_outro1=0;
%         end
%         
%         SM_aqui= SD_aqui + (norm(rDM)/norm(rIQI))*limit_aqui1*gradS_cent_Iele_rIQI;
%         
%         
%         SM_outro= SD_outro + (norm(rDM)/norm(rJQJ))*limit_outro1*gradS_cent_Jele_rJQJ;
% %     end
    
    %% Calculo do fluxo fracional no ponto medio (M)
    
    fe_aqui = flow_fract(visc,SM_aqui,satlimit);
    
    fe_outro = flow_fract(visc,SM_outro,satlimit);
    
    %%   Solver de Reimann pelo Método Upwind
    
    ve_mais = (VN(i) + abs(VN(i)))/2;
    ve_menos = (VN(i) - abs(VN(i)))/2;
    RHS(outro) = RHS(outro) + ve_mais*fe_aqui + ve_menos*fe_outro;
    RHS(aqui)  = RHS(aqui)  - ve_mais*fe_aqui - ve_menos*fe_outro;
end
for iw=1:size(wells,1)
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