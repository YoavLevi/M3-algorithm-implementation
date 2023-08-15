clc
clear all;
close all;
%%
n = 100;
% c = 1/6;
c = 1/65;
max_lambda = 0.9;
max_lambdas = [0.9 0.99];
% max_lambdas = [0.8 0.9 0.95 0.99];
ring_length = uint32(hex2dec('FFFFFFFF'));
iterations = 100;
q_m3_disruption = ceil((n-1).*(max_lambdas./(1-max_lambdas)))+1;
% q_m3_disruption = 2.^ceil(log2(ceil((n-1).*(max_lambdas./(1-max_lambdas)))+1));
q_ACH = zeros(1,iterations);
q_Nginx = zeros(1,iterations);
%% real traces evaluations
tic
% load('hashed_keys315M');
load('unique_keys_hashed_154M');
toc
% hashed_keys = unique(hashed_keys);
% hashed_keys = unique_hashed_keys;
% hashed_keys_sorted = sort(hashed_keys);
toc
req_num = numel(hashed_keys);
real_q_i = zeros(numel(max_lambdas),q_m3_disruption(end));
for k=1:numel(max_lambdas)
    mapping = mod(hashed_keys,q_m3_disruption(k))+1;
    r_i = histcounts(mapping,q_m3_disruption(k));
    real_q_i(k,1:q_m3_disruption(k)) = r_i ./ sum(r_i);
end
%% pre-allocating for speed
max_m3_disruption_lambda_value = zeros(numel(max_lambdas),iterations);
max_m3_disruption_lambda_value_theo = zeros(numel(max_lambdas),iterations);
max_ACH_lambda_value = zeros(1,iterations);
max_ACH_lambda_value_theo = zeros(1,iterations);
max_Ketama_lambda_value = zeros(1,iterations);
max_Ketama_lambda_value_theo = zeros(1,iterations);
%%
for j=1:iterations
    j
    mu_full = randi(10,[1 n]);
%     mu_full = randi(1,[1 n]);
    q_Nginx(j) = sum(mu_full);
    mu = mu_full;
    normalized_mu = mu ./ sum(mu);
    ACH_alloc = ceil(mu./c);
    q = sum(ACH_alloc);
    q_ACH(j) = q;
    while(1)
        vs_keys = char(randi([33 126],q,10));
        vs_keys_Ketama = char(randi([33 126],q_Nginx(j),10));

%         hash_values = crc32_vec(vs_keys);
%         hash_values_Ketama = crc32_vec(vs_keys_Ketama);
        hash_values = fnvhash_vec(vs_keys);
        hash_values_Ketama = fnvhash_vec(vs_keys_Ketama);
        if((numel(hash_values)==numel(unique(hash_values))) && (numel(hash_values_Ketama)==numel(unique(hash_values_Ketama))))
            break;
        end
        fprintf("%g,%g:hash collision, restarting.\n",numel(hash_values),numel(unique(hash_values)));
    end
    
    tic
    for k=1:numel(max_lambdas)
        q_i = return_alloc(q_m3_disruption(k),normalized_mu);
        real_load = return_real_load(q_i,real_q_i(k,1:q_m3_disruption(k)));
        max_m3_disruption_lambda_value(k,j) = inspect_real_stability_func_opt(normalized_mu',real_load');
        max_m3_disruption_lambda_value_theo(k,j) = inspect_stability_func_opt(normalized_mu,q_i);
    end
    toc
    tic
    hash_values_sorted = sort(hash_values);
    hash_values_Ketama_sorted = sort(hash_values_Ketama);
    
    hash_idx = 1;
    hash_idx_Ketama = 1;
    
    count_indices = cumsum(ACH_alloc);
    count_indices_Ketama = cumsum(mu_full);

    scores = zeros(size(mu));
    scores_theo = zeros(size(mu));
    scores_Ketama = zeros(size(mu));
    scores_Ketama_theo = zeros(size(mu));
    
    weights = diff(hash_values_sorted);
    weights = [hash_values_sorted(1)+ring_length-hash_values_sorted(end);weights];

    weights_Ketama = diff(hash_values_Ketama_sorted);
    weights_Ketama = [hash_values_Ketama_sorted(1)+ring_length-hash_values_Ketama_sorted(end);weights_Ketama];

    for i=1:numel(scores)
        while(hash_idx<=count_indices(i))
            idx_on_ring = find(hash_values_sorted==hash_values(hash_idx));
            scores(i) = scores(i) + return_num_req_keys(idx_on_ring,hash_values_sorted,hashed_keys);
%             scores(i) = scores(i) + return_num_req_keys_faster(idx_on_ring,hash_values_sorted,hashed_keys);
            
            % theo
            scores_theo(i) = scores_theo(i) + weights(idx_on_ring);
            
            hash_idx = hash_idx + 1;
        end
        while(hash_idx_Ketama<=count_indices_Ketama(i))
            idx_on_ring_Ketama = find(hash_values_Ketama_sorted==hash_values_Ketama(hash_idx_Ketama));
            scores_Ketama(i) = scores_Ketama(i) + return_num_req_keys(idx_on_ring_Ketama,hash_values_Ketama_sorted,hashed_keys);
%             scores_Ketama(i) = scores_Ketama(i) + return_num_req_keys_opt(idx_on_ring_Ketama,hash_values_Ketama_sorted,hashed_keys);

            % theo
            scores_Ketama_theo(i) = scores_Ketama_theo(i) + weights_Ketama(idx_on_ring_Ketama);

            hash_idx_Ketama = hash_idx_Ketama + 1;
        end
    end

    scores = scores ./ double(req_num);
    scores_Ketama = scores_Ketama ./ double(req_num); 
    
    scores_theo = scores_theo ./ double(ring_length);
    scores_Ketama_theo = scores_Ketama_theo ./ double(ring_length);
    
    max_ACH_lambda_value(j) = inspect_stability_real_ACH_func_opt(scores,normalized_mu);
    max_Ketama_lambda_value(j) = inspect_stability_real_ACH_func_opt(scores_Ketama,normalized_mu);

    max_ACH_lambda_value_theo(j) = inspect_stability_real_ACH_func_opt(scores_theo,normalized_mu);
    max_Ketama_lambda_value_theo(j) = inspect_stability_real_ACH_func_opt(scores_Ketama_theo,normalized_mu);
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
    h(2*k+3) = cdfplot(max_m3_disruption_lambda_value(k,:));
    h(2*k+4) = cdfplot(max_m3_disruption_lambda_value_theo(k,:));
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
% set(h(5),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(3)),')'),'Color','m');
% set(h(6),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(4)),')'),'Color','k');
legend(h([1 3 5 7]),'Location','NorthWest');
% legend(h([1 2 3 4 5 6]),'Location','NorthWest');
title('');
xlim([0 1.005]);
figure(4);
g(1) = cdfplot(q_ACH);
hold on;
g(2) = cdfplot(q_Nginx);
hold on;
for k=1:numel(max_lambdas)
    g(k+2) = cdfplot(q_m3_disruption(k));
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
% set(g(5),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(3)),')'),'LineWidth',1.5,'Color','m');
% set(g(6),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(4)),')'),'LineWidth',1.5,'Color','k');
legend(g([1 2 3 4]),'Location','NorthWest');
% legend(g([1 2 3 4 5 6]),'Location','NorthWest');
% xlim([-inf 105e2]);