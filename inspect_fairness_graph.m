clc
clear all
close all
%% 
q = 10:100;
n=3;
num_sim = 50;
M_mat_int = randi([1 100],num_sim,n);
M_mat = M_mat_int ./ sum(M_mat_int')';
M_mat = M_mat';
lambda = 0.95;
stabilityness = zeros(1,numel(q));
max_load = zeros(1,numel(q));
balance = zeros(1,numel(q));
deviation = zeros(1,numel(q));
for M=M_mat
    for i=1:numel(q)
        [~,~,max_load(i),stabilityness(i),balance(i),deviation(i)] = inspect_fairness_func(q(i),M',lambda);
    end
    figure(5)
    plot(q,max_load);
    hold on;
end
figure(5)
yyaxis left
xlabel('q')
ylabel('max \rho_i')
xline((numel(M)-1) * (lambda/(1-lambda)))
yline(0)
yline(1,'--b')
ylim([lambda 1.15]);
xlim([q(1) inf]);
hold on
alpha_const = lambda*(n-1);
theorem_plot = plot(q,lambda+alpha_const./q,'r','LineWidth',1);
text(q(end)/2,1.005,'Stability Threshold')
legend(theorem_plot,'M3-all');

yyaxis right
axis([-inf inf 1 23/19]);
ylabel('max ^{\rho_i}/_{\rho}');