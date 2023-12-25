clc
clear all
close all
%%
q = 1:1:105;
qmax = q(end);
mu = [0.15 0.23 0.31 0.31]';
max_lambda = zeros(1,numel(q));
for i=1:numel(q)
    max_lambda(i) = inspect_stability_func(mu,q(i));
end
start_index = 1;
figure(2);
plot(q,max_lambda)
xlabel('q')
ylabel('max stable \rho')
hold on
plot(q(start_index:end),(q(start_index:end)./(q(start_index:end)+numel(mu)-1)),'DisplayName','M3 (guaranteed)');
hold on
scatter(q(round(q)==q),max_lambda(round(q)==q),'blue','DisplayName','M3 (measured)');
hold on
grid on;
ax = gca;
ax.XMinorGrid = 'on';
ylim([0.6 1]);
xline(12,':','LineWidth',2,'DisplayName','','Color','#4DBEEE');
xlim([0 105]);
yline(0.8,'--');
legend('','M3 (guaranteed)','M3 (measured)','','','Interpreter','latex','Location','SouthEast','Fontsize',14)
% axes('Position',[0.396214865673831 0.269411691824812 0.349713473088384 0.369828814504302]);
% box on;
% idx = q > 80 & q < 105;
% plot(q(idx),max_lambda(idx));
% hold on;
% plot(q(idx),(q(idx)./(q(idx)+numel(mu)-1)));
% hold on;
% scatter(q(idx),max_lambda(idx),'blue');
% grid on;
% set(gca,'gridlinestyle','--','MinorGridAlpha',0.5)
% axis tight;
% xline(100,'-.','LineWidth',2);
% createrectangles(figure(2));