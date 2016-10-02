function [g]=grad(S,w,esuel,elem,bcflag,bedge,bound)
g=zeros(size(elem,1),3);
for i=1:size(elem,1)
    [m,n]=size(esuel);
    
    for j=1:n
        if esuel(i,j)~=0
            g(i,:)=g(i,:)+ w(:,j,i)'*(S(esuel(i,j))-S(i));
        end
    end
    
end
for jj=1:size(bedge,1);
    g(bedge(jj,3),:)=0;
end
h=1;
b=zeros(100,1);
for ii=1:size(bedge,1);
    
    iele=bedge(ii,3);
    aux=b(:)==iele;
    if max(aux)==0
        t=1;
        list=esuel(iele,:);
        for j=1:length(list)
            if list(j)~=0
                if list(j)==iele
                    bb=bound(iele,t);
                    boundary_type = bedge(bb,5);
                    boundary_condition = bcflag((bcflag(:,1)==boundary_type),2);
                    
                    if boundary_type>200       %neumann
                        
                        p_b =S(iele);
                        
                    elseif boundary_type<200  %dirichlet
                        
                        p_b =2* boundary_condition -S(iele) ;
                    end
                    g(iele,:)=g(iele,:)+ w(:,j,iele)'*(p_b-S(iele));
                    t=t+1;
                    b(h)=iele;
                    
                else
                    
                    g(iele,:)=g(iele,:)+ w(:,j,iele)'*(S(list(j))-S(iele));
                    b(h)=iele;
                end
            end
        end
        
        h=h+1;
    end
end
