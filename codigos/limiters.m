function [phi]=limiters(ielem,ielem1,coord,c,esurn1,esurn2,g,S,elem,volume)
% veja o artigo "Multi-dimensional limiting process for finite volume
% methods on unstructured grids" 2012. pag. 12 e 13
% c: vetor baricentro
% g: matriz gradiente
% S: vetor saturação
% ielem: número de elemento em questão

K1=5; % recomendado pelo author
K2=5; % recomendado pelo author
if elem(ielem,4)~=0 % verifica si o elemento é triangulo o quadrilatero
    j=4;
else
    j=3;
end

for p=1:j
%     delt_menos=dot(g(ielem,:),coord(elem(ielem,p),:)-c(ielem,:));
%     
%    % calcula a saturação max e min en torno a nó em questão
%     [S_max,S_min]=phi_calculate(elem(ielem,p),esurn1,esurn2,S);
%     
%     theta=(S_max-S_min)/(K2*volume(ielem)^1.5); % veja equação 20 do artigo
% 
%     e=(K1/(1+theta))*(S_max-S_min)^2; % veja equaçao 19 do artigo
%     
%     if delt_menos>0 && delt_menos >1e-16
%         delt_mas=S_max-S(ielem);
%         
%         phi(p)=(1/delt_menos)*(((delt_mas^2 +e)*delt_menos +...
%             2*delt_menos^2*delt_mas)/(delt_mas^2 +2*delt_menos^2 +...
%             delt_menos*delt_mas+e));
%         
%     elseif delt_menos<0 && delt_menos>-1e-16
%         delt_mas= S_min-S(ielem);
%                 
%         phi(p)=(1/delt_menos)*(((delt_mas^2 +e)*delt_menos +...
%             2*delt_menos^2*delt_mas)/(delt_mas^2 +2*delt_menos^2 +...
%             delt_menos*delt_mas+e));
%         
%     else
%         phi(p)=1;
%     end


        delt_menos=dot(g(ielem,:),coord(elem(ielem,p),:)-c(ielem,:));
        [S_max,S_min]=phi_calculate(elem(ielem,p),esurn1,esurn2,S);
        
        theta=(S_max-S_min)/(K2*volume(ielem)^1.5);
        
        e=(K1/(1+theta))*(S_max-S_min)^2;
        
        if (delt_menos>0 && delt_menos >1e-16)|| (delt_menos<0 && delt_menos>-1e-16)
            if ((S_max-S(ielem))/delt_menos)>((S_min-S(ielem))/delt_menos)
                delt_mas= S_max-S(ielem);
            else
                delt_mas= S_min-S(ielem);
            end
            
            phi(p)=(1/delt_menos)*(((delt_mas^2 +e)*delt_menos + 2*delt_menos^2*delt_mas)/(delt_mas^2 +2*delt_menos^2 + delt_menos*delt_mas+e));
           
        else
            phi(p)=1;
        end
  
end
end