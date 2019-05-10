

function Int=Lyapunov(e1,e2,phi,x)
global P1 P2 P3 mu
Y=zeros(1,1);
for i=1:max(size(e1))
    
   for j=1:min(size(e1))  
    Pz=[P1*exp(-mu*x(j)),P2'; P2, P3];
   Y(j)=[e1(i,j), e2(i,j), phi(i,:)]*Pz*[e1(i,j), e2(i,j),  phi(i,:)]';
    end
    Int(i)=(trapz(x,Y));
end
end