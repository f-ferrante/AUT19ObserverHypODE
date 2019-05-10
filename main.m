%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Francesco Ferrante, Gipsa-lab & Université Grenoble Alpes 
% Yeat: 2019
%
% Description: This code enables to reproduce the example in the paper by
% F., Cristofaro, Prieur, 2019. 

% This code requires YALMIP parser for Linear Matrix Inequality, freely avaialbe at https://yalmip.github.io. 
% Any SDP solver can be used. Here we used SDPT3 freely avaialbe at https://github.com/SQLP/SDPT3

% The paremeters mu and theta needs to be suitably tuned. 

% The code determines the observer gain and simulates the response of the
% interconnection of the plant and the observer from the initial condition given in the paper. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;

Thorizon=100;  %Simulation horizon

%System paremeters
Lambda=[sqrt(2),0;0,2];
F=0.5*[-1,0.2;1,0.2];
M=[1, 1];
A=[0 1 0; 0 0 1; 0 0 0];
C=[1 0 0;0 0 1];
n_x=max(size(Lambda));
n_chi=max(size(A));
n_y=min(size(M));

%Desired convergence rate
lc=0.03;
%Select values for mu and vartheta
vartheta=1;
mu=1.3;
%Define SDP problem 
P1=diag(sdpvar(n_x,1));
P2=sdpvar(n_chi,n_x,'full');
P3=sdpvar(n_chi,n_chi,'symmetric');
T=sdpvar(n_chi,n_chi,'symmetric');
J=sdpvar(n_chi,n_y,'full');

Upsilon=[exp(-mu)*(-mu*Lambda*P1-F'*P1-P1*F), zeros(n_x, n_x), P2'*A-F'*P2';
         zeros(n_x, n_x), -Lambda*P1*exp(-mu), -Lambda*P2'-M'*J';
   (P2'*A-F'*P2')', (-Lambda*P2'-M'*J')', (P3*A+P2*Lambda*C)+(P3*A+P2*Lambda*C)'+C'*Lambda*P1*C]...
         +lc*[P1*exp(-mu), zeros(n_x, n_x), P2'; 
              zeros(n_x,n_x), zeros(n_x,n_x), zeros(n_x,n_chi); 
              P2,zeros(n_x,n_chi)', P3];

Gamma=[-P2'; zeros(n_x,n_chi); zeros(n_chi,n_chi)];

Omega=[zeros(n_y,n_x),M,zeros(n_y,n_chi)];

Sigma=[Upsilon, Gamma, Omega'*J'; 
    Gamma', -vartheta*P3, zeros(n_chi, n_chi); 
    (Omega'*J')',zeros(n_chi, n_chi)', -1/vartheta*P3];

P=[P1*exp(-mu), P2'; P2, P3];

%Solve SDP problem
problem=[Sigma<=0, P>=1e-6*eye(n_chi+n_x), P<=3000*eye(n_chi+n_x)];
options=sdpsettings('solver','sdpt3','verbose',3);
solution=solvesdp(problem,0,options);
%Define paremeter for simulation
if(solution.problem==0)
warning off;
L=inv(double(P3))*double(J);
P1=double(P1);
P2=double(P2);
P3=double(P3);
Fr=F;
lamb1=Lambda(1,1);
lamb2=Lambda(2,2);
global L P1 P2 P3 A Fr mu C M lamb1 lamb2 lc
%Run Simulation
runSim(Thorizon);
else
    disp('It looks like your problem is infeasible');
end

