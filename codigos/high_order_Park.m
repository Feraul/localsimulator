function [S,limiter_lef,limiter_rel] = high_order_Park(S,VN,dt,volume,satlimit,inedge,...
    c,visc,elem,xM,bedge,coord,pormap,wells,q,f_elem,bcflag,esurn1,esurn2,Sat_max_min,w,bound,FF,esuel)

% c: baricentro da célula
% gradient in element
%g = grad_S(elem,c,S,esuel1,esuel2,local_volume,viz_b,bedge,bcflag,inedge);
g=grad(S,w,esuel,elem,bcflag,bedge,bound);
RHS = zeros(size(S,1),1);
limiter_rel=zeros(size(inedge,1),1);
limiter_lef=zeros(size(inedge,1),1);
for i = 1:size(inedge,1)
    
    v1 = inedge(i,1);
    v2 = inedge(i,2);
    
    rel=inedge(i,4);    %indice do elemento a direita da face i
    lef=inedge(i,3);    %indice do elemento a esquerda da face i
    por=pormap(1);
    
    %% limiter 1
    [phi_lef]=limiter_woodfield(Sat_max_min,lef,S);
    [phi_rel]=limiter_woodfield(Sat_max_min,rel,S);
    %% limiter 2
    [tau_lef]=limiters(lef,rel,coord,c,esurn1,esurn2,g,S,elem,volume);
    
    [tau_rel]=limiters(rel,lef,coord,c,esurn1,esurn2,g,S,elem,volume);
    
%     switch HO
%         case {'Park','Venkatakrishnan'}
%             
            
            limiter_lef(i)=min(tau_lef);
            limiter_rel(i)=min(tau_rel);
            
            SM_L = S(lef) + limiter_lef(i)*dot(g(lef,:), xM(i,:)-c(lef,:));
            SM_R = S(rel) + limiter_rel(i)*dot(g(rel,:), xM(i,:)-c(rel,:));
            
%         case 'Hirsch'
%             limiter_lef(i)= tau_lef;
%             limiter_rel(i)=tau_rel;
%             SM_L = S(lef) + phi_lef*limit1(i)*0.25*((1-k)*gradS_upwind_aqui+(1+k)*gradS_cent_aqui);
%             SM_R = S(rel) + phi_rel*limit2(i)*0.25*((1-k)*gradS_upwind_outro+(1+k)*gradS_cent_outro);
%     end
    
    [fe_L] = flow_fract(visc,SM_L,satlimit); % f(SM_L)
    
    [fe_R] = flow_fract(visc,SM_R,satlimit); % f(SM_R)
    
    
    ve_mais = (VN(i) + abs(VN(i)))/2;
    ve_menos = (VN(i) - abs(VN(i)))/2;
    RHS(rel) = RHS(rel) + ve_mais*fe_L + ve_menos*fe_R;
    RHS(lef)  = RHS(lef)  - ve_mais*fe_L - ve_menos*fe_R;
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