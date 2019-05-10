function runSim(Thorizon)
global lc;
[zg,tg, x1, x2, xhat1, xhat2, chi1, chi2, chi3, chihat1, chihat2, chihat3]=hyp_dynamic_bc(Thorizon);
z=zg(1,:);
eta=[chi1-chihat1,chi2-chihat2,chi3-chihat3];
W=Lyapunov(x1-xhat1,x1-xhat1,eta,z);
hold on
grid on
figure(1)
plot(tg, W./max(W),'-b','linewidth', 2);
xlabel('$t$','Interpreter','latex');
title('Normalized Lyapunov Functional');
plot(tg, exp(-lc*tg), ':k','linewidth', 2); 
%%%
figure(2)
hold on
grid on
colormap(winter);
hold on
subplot(2,1,1)
mesh(tg,zg,x1-xhat1);
title('$\varepsilon_1$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
subplot(2,1,2)
mesh(tg,zg,x2-xhat2);
title('$\varepsilon_2$','Interpreter','latex');
xlabel('$t$','Interpreter','latex');
ylabel('$z$','Interpreter','latex');
%%%
figure(3)
hold on
grid on
subplot(3,1,1)
plot(tg, chi1-chihat1,'-b','linewidth', 2)
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$\eta_1$','Interpreter','latex');
subplot(3,1,2)
plot(tg, chi2-chihat2,'-b','linewidth', 2)
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$\eta_2$','Interpreter','latex');
subplot(3,1,3)
plot(tg, chi3-chihat3,'-b','linewidth', 2)
grid on
xlabel('$t$','Interpreter','latex');
ylabel('$\eta_3$','Interpreter','latex');
end
