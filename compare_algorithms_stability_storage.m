clc
clear all;
close all;
%%
n = 30;
c = 1;
max_lambda = 0.9;
max_lambdas = [0.9 0.99];
ring_length = uint32(hex2dec('FFFFFFFF'));
iterations = 1000;
q_m3 = zeros(1,iterations);
q_m3_disruption = floor((n-1).*(max_lambdas./(1-max_lambdas)))+1;
q_ACH = zeros(1,iterations);
q_Nginx = zeros(1,iterations);
servers_weights = [2 5];
%% real traces evaluations
filename = 'testDB.txt';
lines = readlines(filename);
values = split(lines(1:end-1),',');
r = char(join(values(2:end,3:7),''));
hashed_keys = crc32_vec(r);
req_num = numel(hashed_keys);
for k=1:numel(max_lambdas)
    mapping = mod(hashed_keys,q_m3_disruption(k))+1;
    r_i = histcounts(mapping,q_m3_disruption(k));
    if k==1
        real_q_i_first = r_i ./ sum(r_i);
    elseif k==2
        real_q_i_second = r_i ./ sum(r_i);
    end
end
%% pre-allocating for speed
max_m3_disruption_lambda_value = zeros(numel(max_lambdas),iterations);
max_ACH_lambda_value = zeros(1,iterations);
max_Ketama_lambda_value = zeros(1,iterations);
%%
for j=1:iterations
    j
    heavy_count = randi(15);
    small_count = randi(15);
    mu_full = [ones(1,heavy_count).*servers_weights(2),ones(1,small_count).*servers_weights(1)];
    q_Nginx(j) = sum(mu_full);
    mu = 100 * mu_full;
    normalized_mu = mu ./ sum(mu);
    ACH_alloc = ceil(mu./c);
    q = sum(ACH_alloc);
    q_ACH(j) = q;
%     hash_values = zeros(1,q);
    while(1)
        vs_keys = char(randi([33 126],q,10));
        vs_keys_Ketama = char(randi([33 126],q_Nginx(j),10));

        hash_values = crc32_vec(vs_keys);
        hash_values_Ketama = crc32_vec(vs_keys_Ketama);
        if((numel(hash_values)==numel(unique(hash_values))) && (numel(hash_values_Ketama)==numel(unique(hash_values_Ketama))))
            break;
        end
        fprintf("%g,%g:hash collision, restarting.\n",numel(hash_values),numel(unique(hash_values)));
    end
    
    for k=1:numel(max_lambdas)
        if k==1
            real_q_i = real_q_i_first;
        elseif k==2
            real_q_i = real_q_i_second;
        end
        q_i = return_alloc(q_m3_disruption(k),normalized_mu);
        real_load = return_real_load(q_i,real_q_i);
        max_m3_disruption_lambda_value(k,j) = inspect_real_stability_func(normalized_mu',real_load');
    end
    hash_values_sorted = sort(hash_values);
    hash_values_Ketama_sorted = sort(hash_values_Ketama);
    
    hash_idx = 1;
    hash_idx_Ketama = 1;
    
    count_indices = cumsum(ACH_alloc);
    count_indices_Ketama = cumsum(ACH_alloc./(q/q_Nginx(j)));

    scores = zeros(size(mu));
    scores_Ketama = zeros(size(mu));
  
    for i=1:numel(scores)
        while(hash_idx<=count_indices(i))
            idx_on_ring = find(hash_values_sorted==hash_values(hash_idx));
            scores(i) = scores(i) + return_num_req_keys(idx_on_ring,hash_values_sorted,hashed_keys);
            hash_idx = hash_idx + 1;
        end
        while(hash_idx_Ketama<=count_indices_Ketama(i))
            idx_on_ring_Ketama = find(hash_values_Ketama_sorted==hash_values_Ketama(hash_idx_Ketama));
            scores_Ketama(i) = scores_Ketama(i) + return_num_req_keys(idx_on_ring_Ketama,hash_values_Ketama_sorted,hashed_keys);
            hash_idx_Ketama = hash_idx_Ketama + 1;
        end
    end
    scores = scores ./ double(req_num);
    scores_Ketama = scores_Ketama ./ double(req_num); 
    
    max_ACH_lambda_value(j) = inspect_stability_real_ACH_func(scores,normalized_mu);
    max_Ketama_lambda_value(j) = inspect_stability_real_ACH_func(scores_Ketama,normalized_mu);
end
markers_res = 20;
figure(1);
h(3) = cdfplot(max_Ketama_lambda_value);
hold on
h(1) = cdfplot(max_ACH_lambda_value);
hold on
for k=1:numel(max_lambdas)
    h(k+3) = cdfplot(max_m3_disruption_lambda_value(k,:));
end
xlabel('max stable \rho');
ylabel('CDF');
set( h(1), 'LineStyle', '--', 'Color', '#EDB120','DisplayName','ACH','Marker','o','MarkerIndices',1:floor(numel(h(1).XData)/markers_res):numel(h(1).XData/50));
set( h(3), 'Color', 'red','DisplayName','NGINX','Marker','*','MarkerIndices',1:floor(numel(h(3).XData)/markers_res):numel(h(3).XData));
set(h(4),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(1)),')'),'LineWidth',1.5,'Color','#0072BD','Marker','+','MarkerIndices',1:floor(numel(h(4).XData)/markers_res):numel(h(4).XData));
set(h(5),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(2)),')'),'LineWidth',1.5,'Color','#77AC30');
legend(h([1 3 4 5]),'Location','NorthWest');
title('');
xlim([0.2 inf]);
figure(4);
g(1) = cdfplot(q_ACH);
hold on;
hold on;
g(3) = cdfplot(q_Nginx);
hold on;
for k=1:numel(max_lambdas)
    g(k+3) = cdfplot(q_m3_disruption(k));
end
hold on;
set( g(1), 'Color', '#EDB120','DisplayName','ACH','LineStyle','--','Marker','o','MarkerIndices',1:floor(numel(g(1).XData)/markers_res):numel(g(1).XData));
set( g(3), 'Color', 'red','DisplayName','NGINX','Marker','*','MarkerIndices',1:floor(numel(g(3).XData)/markers_res):numel(g(3).XData));
title('');
set(gca,'XScale','log')
xlabel('q');
ylabel('CDF');

set(g(4),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(1)),')'),'LineWidth',1.5,'Color','#0072BD','Marker','+');
set(g(5),'DisplayName',strcat('M3-all (max \rho=',num2str(max_lambdas(2)),')'),'LineWidth',1.5,'Color','#77AC30');
