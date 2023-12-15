clc
clear all;
close all;
%%
n = 100;
c = 1/65;
max_lambdas = [0.9 0.99];
ring_length = uint32(hex2dec('FFFFFFFF'));
num_of_iterations = 100;
q_m3 = ceil((n-1).*(max_lambdas./(1-max_lambdas)))+1;
q_ACH = zeros(1,num_of_iterations);
q_Nginx = zeros(1,num_of_iterations);
%% real traces evaluations
tic
load('unique_keys_hashed_154M');
toc
req_num = numel(hashed_keys);
real_q_i = zeros(numel(max_lambdas),q_m3(end));
for k=1:numel(max_lambdas)
    mapping = mod(hashed_keys,q_m3(k))+1;
    r_i = histcounts(mapping,q_m3(k));
    real_q_i(k,1:q_m3(k)) = r_i ./ sum(r_i);
end
%% pre-allocating for speed
max_m3_lambda_value = zeros(numel(max_lambdas),num_of_iterations);
max_m3_lambda_value_theo = zeros(numel(max_lambdas),num_of_iterations);
max_ACH_lambda_value = zeros(1,num_of_iterations);
max_ACH_lambda_value_theo = zeros(1,num_of_iterations);
max_Ketama_lambda_value = zeros(1,num_of_iterations);
max_Ketama_lambda_value_theo = zeros(1,num_of_iterations);
%% Main loop for each simulation instance
for iter=1:num_of_iterations
    fprintf("simulation No. %g\n",iter);
    mu_full = randi(10,[1 n]);
    q_Nginx(iter) = sum(mu_full);
    mu = mu_full;
    normalized_mu = mu ./ sum(mu);
    ACH_alloc = ceil(mu./c);
    q = sum(ACH_alloc);
    q_ACH(iter) = q;
    while(1)
        vs_keys = char(randi([33 126],q,10));
        vs_keys_Ketama = char(randi([33 126],q_Nginx(iter),10));
        hash_values = fnvhash_vec(vs_keys);
        hash_values_Ketama = fnvhash_vec(vs_keys_Ketama);
        if((numel(hash_values)==numel(unique(hash_values))) && (numel(hash_values_Ketama)==numel(unique(hash_values_Ketama))))
            break;
        end
        fprintf("%g,%g:hash collision, restarting.\n",numel(hash_values),numel(unique(hash_values))); % hitting this line is expected
    end
    for k=1:numel(max_lambdas)
        q_i = return_alloc(q_m3(k),normalized_mu);
        real_load = return_real_load(q_i,real_q_i(k,1:q_m3(k)));
        max_m3_lambda_value(k,iter) = inspect_real_stability_func_opt(normalized_mu',real_load');
        max_m3_lambda_value_theo(k,iter) = inspect_stability_func_opt(normalized_mu,q_i);
    end
    tic
    hash_values_sorted = sort(hash_values);
    hash_values_Ketama_sorted = sort(hash_values_Ketama);
    
    hash_idx_ACH = 1;
    hash_idx_Ketama = 1;
    
    count_indices = cumsum(ACH_alloc);
    count_indices_Ketama = cumsum(mu_full);

    scores_ACH = zeros(size(mu));
    scores_ACH_theo = zeros(size(mu));
    scores_Ketama = zeros(size(mu));
    scores_Ketama_theo = zeros(size(mu));
    
    weights_ACH = diff(hash_values_sorted);
    weights_ACH = [hash_values_sorted(1)+ring_length-hash_values_sorted(end);weights_ACH];

    weights_Ketama = diff(hash_values_Ketama_sorted);
    weights_Ketama = [hash_values_Ketama_sorted(1)+ring_length-hash_values_Ketama_sorted(end);weights_Ketama];

    for i=1:numel(scores_ACH)
        while(hash_idx_ACH<=count_indices(i))
            idx_on_ring = find(hash_values_sorted==hash_values(hash_idx_ACH));
            scores_ACH(i) = scores_ACH(i) + return_num_req_keys(idx_on_ring,hash_values_sorted,hashed_keys);
            scores_ACH_theo(i) = scores_ACH_theo(i) + weights_ACH(idx_on_ring);
            hash_idx_ACH = hash_idx_ACH + 1;
        end
        while(hash_idx_Ketama<=count_indices_Ketama(i))
            idx_on_ring_Ketama = find(hash_values_Ketama_sorted==hash_values_Ketama(hash_idx_Ketama));
            scores_Ketama(i) = scores_Ketama(i) + return_num_req_keys(idx_on_ring_Ketama,hash_values_Ketama_sorted,hashed_keys);
            scores_Ketama_theo(i) = scores_Ketama_theo(i) + weights_Ketama(idx_on_ring_Ketama);
            hash_idx_Ketama = hash_idx_Ketama + 1;
        end
    end
    scores_ACH = scores_ACH ./ double(req_num);
    scores_Ketama = scores_Ketama ./ double(req_num); 
    scores_ACH_theo = scores_ACH_theo ./ double(ring_length);
    scores_Ketama_theo = scores_Ketama_theo ./ double(ring_length);
    max_ACH_lambda_value(iter) = inspect_stability_real_ACH_func_opt(scores_ACH,normalized_mu);
    max_Ketama_lambda_value(iter) = inspect_stability_real_ACH_func_opt(scores_Ketama,normalized_mu);
    max_ACH_lambda_value_theo(iter) = inspect_stability_real_ACH_func_opt(scores_ACH_theo,normalized_mu);
    max_Ketama_lambda_value_theo(iter) = inspect_stability_real_ACH_func_opt(scores_Ketama_theo,normalized_mu);
    toc
end
markers_res = 20;
figure(1);
h(1) = cdfplot(max_ACH_lambda_value);
hold on
h(2) = cdfplot(max_ACH_lambda_value_theo);
hold on
h(3) = cdfplot(max_Ketama_lambda_value);
hold on
h(4) = cdfplot(max_Ketama_lambda_value_theo);
hold on
for k=1:numel(max_lambdas)
    h(2*k+3) = cdfplot(max_m3_lambda_value(k,:));
    h(2*k+4) = cdfplot(max_m3_lambda_value_theo(k,:));
end
xlabel('max stable \rho');
ylabel('CDF');
set(h(1), 'Color', '#EDB120','DisplayName','ACH','Marker','o','MarkerIndices',1:floor(numel(h(1).XData)/markers_res):numel(h(1).XData/50));
set(h(2), 'LineStyle', '--', 'Color', '#EDB120','DisplayName','ACH','Marker','o','MarkerIndices',1:floor(numel(h(2).XData)/markers_res):numel(h(2).XData/50));
set(h(3), 'Color', 'red','DisplayName','NGINX','Marker','*','MarkerIndices',1:floor(numel(h(3).XData)/markers_res):numel(h(3).XData));
set(h(4), 'LineStyle', '--','Color', 'red','DisplayName','NGINX','Marker','*','MarkerIndices',1:floor(numel(h(4).XData)/markers_res):numel(h(4).XData));
set(h(5), 'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(1)),')'),'LineWidth',1.5,'Color','#0072BD','Marker','+','MarkerIndices',1:floor(numel(h(5).XData)/markers_res):numel(h(5).XData));
set(h(6), 'LineStyle', '--','DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(1)),')'),'LineWidth',1.5,'Color','#0072BD','Marker','+','MarkerIndices',1:floor(numel(h(6).XData)/markers_res):numel(h(6).XData));
set(h(7), 'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(2)),')'),'LineWidth',1.5,'Color','#77AC30');
set(h(8), 'LineStyle', '--','DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(2)),')'),'LineWidth',1.5,'Color','#77AC30');
legend(h([1 3 5 7]),'Location','NorthWest');
title('');
xlim([0 1.005]);
figure(4);
g(1) = cdfplot(q_ACH);
hold on;
g(2) = cdfplot(q_Nginx);
hold on;
for k=1:numel(max_lambdas)
    g(k+2) = cdfplot(q_m3(k));
end
hold on;
set( g(1), 'Color', '#EDB120','DisplayName','ACH','Marker','o','MarkerIndices',1:floor(numel(g(1).XData)/markers_res):numel(g(1).XData));
set( g(2), 'Color', 'red','DisplayName','NGINX','Marker','*','MarkerIndices',1:floor(numel(g(3).XData)/markers_res):numel(g(3).XData));
title('');
set(gca,'XScale','log')
xlabel('q');
ylabel('CDF');
set(g(3),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(1)),')'),'LineWidth',1.5,'Color','#0072BD','Marker','+');
set(g(4),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(2)),')'),'LineWidth',1.5,'Color','#77AC30');
legend(g([1 2 3 4]),'Location','NorthWest');