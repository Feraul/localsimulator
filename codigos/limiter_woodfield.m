function [phi]=limiter_woodfield(Sat_max_min,ielem,S)

if Sat_max_min(ielem,1)-Sat_max_min(ielem,2)<=1e-20
    phi= 1;
else
    gamma= (S(ielem)-Sat_max_min(ielem,2))/(Sat_max_min(ielem,1)-Sat_max_min(ielem,2));
    if gamma >=1 | gamma<=0
        phi=0;
    elseif 0.2<= gamma && gamma<=1-0.2
        phi=1;
    elseif 0< gamma && gamma< 0.2
        phi=gamma/0.2;
    elseif 1-0.2 < gamma && gamma<1
        phi=(1-gamma)/0.2;
    end
end
end