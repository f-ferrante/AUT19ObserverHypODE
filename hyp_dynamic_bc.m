
function [x1,t1, V1to, V2to,V3to, V4to, chi1to,chi2to,chi3to, chi1_hatto, chi2_hatto, chi3_hatto]=hyp_dynamic_bc(Thorizon)
fprintf('Completed: %3.0f%%',0);
global lamb1 lamb2 A C Fr L M  

method = 'LxF';
%method = 'SLxW';
pdefun = @pdes;

% Mesh grid
Nx = 10; % number of points in the space domain   
x = linspace(0,1,Nx);
Tinit = 0;
t = Tinit;

% Initial conditions for x dynamics
V(1,:)=0.1*(cos(2*pi*x)-1);
V(2,:)=-0.1*(cos(4*pi*x)-1);
% Initial conditions for \hat{x} dynamics
V(3,:)=0;
V(4,:)=0;
% intial condition for the chi dynamics
chi10=0.1;
chi20=-0.1;
chi30=0.2;
V(5,:)=chi10;
V(6,:)=chi20;
V(7,:)=chi30;
% intial condition for the \hat{chi} dynamics
chi1_hat0=0;
chi2_hat0=0;
chi3_hat0=0;
V(8,:)=chi1_hat0;
V(9,:)=chi2_hat0;
V(10,:)=chi3_hat0;
sol = setup(1,pdefun,t,x,V,method,[],@bcfun);
% check CFL condition
dx = x(2)-x(1);
dt = 0.95*dx/20;
% Simulation completion indeicators
Nt=floor(Thorizon/dt);
howfar = Thorizon/Nt;
% to stack the solution
V1to =zeros(Nt,Nx);
V2to =zeros(Nt,Nx);
V3to =zeros(Nt,Nx);
V4to =zeros(Nt,Nx);
chi1to =zeros(Nt,1);
chi2to =zeros(Nt,1);
chi3to =zeros(Nt,1);
chi1_hatto =zeros(Nt,1);
chi2_hatto =zeros(Nt,1);
chi3_hatto =zeros(Nt,1);
%Intial conditions are fed into the simulator
chi1to(1,1)=chi10;
chi2to(1,1)=chi20;
chi3to(1,1)=chi30;
chi1_hatto(1,1)=chi1_hat0;
chi2_hatto(1,1)=chi2_hat0;
chi3_hatto(1,1)=chi3_hat0;
V1to(1,:) = V(1,:);
V2to(1,:) = V(2,:);
V3to(1,:) = V(3,:);
V4to(1,:) = V(4,:);
chi1to(1,:) = V(5,1);
chi2to(1,:) = V(6,1);
chi3to(1,:) = V(7,1);
chi1_hatto(1,:) = V(8,1);
chi2_hatto(1,:) = V(9,1);
chi3_hatto(1,:) = V(10,1);
%Solution to the PDEs
for m = 1:Nt-1
    indicator=floor(m/(Nt-1)*100);
   % ['Simulation: ' num2str(m) '/' num2str(Nt)]
    %disp(['Simulation:', num2str(indicator), '%']);
    fprintf('\b\b\b\b%3.0f%%',indicator);
    sol = hpde(sol,howfar,dt);  
    t = sol.t;
    V = sol.u;
    V1to(m+1,:) = V(1,:);
    V2to(m+1,:) = V(2,:);
    V3to(m+1,:) = V(3,:);
    V4to(m+1,:) = V(4,:);
    chi1to(m+1) = V(5,1);
    chi2to(m+1) = V(6,1);
    chi3to(m+1) = V(7,1);
    chi1_hatto(m+1) = V(8,1);
    chi2_hatto(m+1) = V(9,1);
    chi3_hatto(m+1) = V(10,1);
end
t = linspace(Tinit,Thorizon,Nt);
[x1,t1] = meshgrid(x,t);
end


    
%=========================================================================
% Subfunctions
    
    function F = pdes(t,x,V,V_x) % linear PDE
         global lamb1 lamb2 A Fr L M C
         %keyboard 
          F = zeros(size(V));
  
        F(1,:) = -lamb1.*V_x(1,:)-Fr(1,1).*V(1,:)-Fr(1,2).*V(2,:);
        F(2,:) = -lamb2.*V_x(2,:)-Fr(2,1).*V(1,:)-Fr(2,2).*V(2,:);
        
        F(3,:) = -lamb1.*V_x(3,:)-Fr(1,1).*V(3,:)-Fr(1,2).*V(4,:);
        F(4,:) = -lamb2.*V_x(4,:)-Fr(2,1).*V(3,:)-Fr(2,2).*V(4,:);
        
        chidot=A*[V(5,1);V(6,1);V(7,1)];
        
        chidot_hat=A*[V(8,1);V(9,1);V(10,1)]+L*(M*[V(1,end);V(2,end)]-M*[V(3,end);V(4,end)]);
      
        F(5,:) = chidot(1);
        F(6,:) = chidot(2);
        F(7,:) = chidot(3);
        
        F(8,:) = chidot_hat(1);
        F(9,:) = chidot_hat(2);
        F(10,:) = chidot_hat(3);
    end
  
        % end function pdes
    
    function [XL,XR] = bcfun(t,XLex,XRex)
    global C
    XLtmp = XLex;
    % dynamic boundary condition with eta
        %keyboard
       % Xbc=0*[XRex(1);XRex(2)]+C*[XRex(3);XRex(4)];
        Xbc=C*[XRex(5);XRex(6);XRex(7)];
        Xbc_hat=C*[XRex(8);XRex(9);XRex(10)];
        %disp(['boundary condition for X_1 and X_2:' num2str(Xbc')])
        XLtmp(1) = Xbc(1);
        XLtmp(2) = Xbc(2);
        XLtmp(3) = Xbc_hat(1);
        XLtmp(4) = Xbc_hat(2);
        XRtmp = XRex;
    XL=XLtmp;
    XR=XRtmp;
    end
