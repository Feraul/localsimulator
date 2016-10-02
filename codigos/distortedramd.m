function [auxcoord]=distortedramd
global coord bedge

auxcoord=zeros(size(coord,1),4);

a=-0.5;
b=0.5;

ex = a + (b-a).*rand(size(coord,1),1);
ey = a + (b-a).*rand(size(coord,1),1);
alpha=0.55; 
h=1/192;

for icoord=1:size(coord,1)
   
    if icoord>size(bedge,1) && abs(coord(icoord,1)-0.5)>1e-10
        x=coord(icoord,1)+alpha*ex(icoord,1)*h;
        
        y=coord(icoord,2)+alpha*ey(icoord,1)*h;
        auxcoord(icoord,1)=icoord;
        auxcoord(icoord,2)=x;
        auxcoord(icoord,3)=y;
    else
        x=coord(icoord,1);
        
        y=coord(icoord,2);
        auxcoord(icoord,1)=icoord;
        auxcoord(icoord,2)=x;
        auxcoord(icoord,3)=y;
    end
end

end
