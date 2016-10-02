function [auxphi]=limiterwoodf(Sat_max_min, ielem,upsilon, S_old)

if Sat_max_min(ielem,1)-Sat_max_min(ielem,2)<1e-20
    auxphi= 1;
else
    gamma= (S_old(ielem)-Sat_max_min(ielem,2))/(Sat_max_min(ielem,1)-Sat_max_min(ielem,2));
    if (1<gamma || gamma==1) || (gamma<0 || gamma==0)
        auxphi=0;
    elseif (upsilon<gamma ||gamma==upsilon) && (gamma<1-upsilon || gamma==1-upsilon)
        auxphi=1;
    elseif 0< gamma && gamma<upsilon
        auxphi=gamma/upsilon;
    elseif 1-upsilon < gamma && gamma<1
        auxphi=(1-gamma)/upsilon;
    end
end

end