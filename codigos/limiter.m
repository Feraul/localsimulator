function L_VA = limiter(a,b,volume)

%Van Albada limiting function

 e = 1e-16;
%e=sqrt((volume*0.3)^3);
%beta= 1; % quado beta=1 ---> MinMod (o patamar é um poqueno menor que de Van Albada) e quando beta=2 ---> SuperBee (o patamar ainda maior)
 csi = (b^2 + e^2)/(a^2 + b^2 + 2*e^2);
 
if (a*b)>0
    L_VA = csi*a + (1-csi)*b;
      
else
    L_VA = 0;
end 
%L_VA=max(0, (2*a*b+e)/(a^2 + b^2 +e)); % limitador lonher
%L_VA= b*max([0,min(beta*a/b,1),min(a/b,beta)]);
 
end