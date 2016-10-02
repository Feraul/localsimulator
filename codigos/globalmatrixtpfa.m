function [M,I] = globalmatrixtpfa(Kde, Kn, nflag, Hesq,wells,mobility)

global inedge bedge elem coord centelem bcflag

% Constrói a matriz global.
% prealocação da matriz global e do vetor termo de fonte
M=zeros(size(elem,1),size(elem,1));
I=zeros(size(elem,1),1);
auxmobility1=mobility(1:size(inedge,1),1);
auxmobility2=mobility((size(inedge,1)+1):(size(inedge,1)+size(bedge,1)),1);
mobility(1:size(bedge,1),1)=auxmobility2;
mobility((size(bedge,1)+1):(size(inedge,1)+size(bedge,1)),1)=auxmobility1;

if max(wells)~=0
    sumvol=0;
    for iw = 1:size(wells,1)
        
        if wells(iw,2)==1            % injetor
            I(wells(iw,1))= 1*centelem(wells(iw,1));        % injeta um m3 de agua por dia (d)
            sumvol=sumvol+ centelem(wells(iw,1));
        end
    end
    I=I./sumvol;
else
    
    % Loop de faces de contorno
    
    for ifacont=1:size(bedge,1)
        v0=coord(bedge(ifacont,2),:)-coord(bedge(ifacont,1),:);
        
        % calculo das constantes nas faces internas
        A=-Kn(ifacont)/(Hesq(ifacont)*norm(v0));     
        
        if bedge(ifacont,5)<200
            
            c1=nflag(bedge(ifacont,1),2);
            c2=nflag(bedge(ifacont,2),2);
            
            %Preenchimento
            
            M(bedge(ifacont,3),bedge(ifacont,3))=M(bedge(ifacont,3),bedge(ifacont,3))- mobility(ifacont)*A*(norm(v0)^2);
            
            I(bedge(ifacont,3))=I(bedge(ifacont,3))-mobility(ifacont)*c1*A*(norm(v0)^2);
            
        else
            x=bcflag(:,1)==bedge(ifacont,5);
            r=find(x==1);
              I(bedge(ifacont,3))=I(bedge(ifacont,3))- normcont*bcflag(r,2);
          
            
        end
    end
end
for iface=1:size(inedge,1),
    
    %Contabiliza as contribuições do fluxo numa faces  para os elementos %
    %a direita e a esquerda dela.                                        %
    M(inedge(iface,3), inedge(iface,3))=M(inedge(iface,3), inedge(iface,3))-mobility(iface+size(bedge,1))*Kde(iface,1);
    M(inedge(iface,3), inedge(iface,4))=M(inedge(iface,3), inedge(iface,4))+mobility(iface+size(bedge,1))*Kde(iface,1);
    M(inedge(iface,4), inedge(iface,4))=M(inedge(iface,4), inedge(iface,4))-mobility(iface+size(bedge,1))*Kde(iface,1);
    M(inedge(iface,4), inedge(iface,3))=M(inedge(iface,4), inedge(iface,3))+mobility(iface+size(bedge,1))*Kde(iface,1);
end
% Fim do loop de faces internas
for iw = 1:size(wells,1)
    if wells(iw,2)==2 %produtor
        M(wells(iw,1),:)=0*M(wells(iw,1),:);
        M(wells(iw,1),wells(iw,1))=1;
        I(wells(iw,1))=0;
    end
end
end
