function [max_lambda_value] = inspect_stability_real_ACH_func_opt(rl_ACH,normalized_mu)
%     res = 1e-5;
%     lambda = 0:res:1-res;
    lambda = 1;
%     rl_ACH = alloc' .* lambda ./ (q.*normalized_mu)';
    valid_indices = find(max((rl_ACH'.*lambda)./(normalized_mu')) <= 1);
    rho_max = max((rl_ACH'.*lambda)./(normalized_mu'));
%     max_lambda_value = lambda(valid_indices(end));
    max_lambda_value = 1/rho_max;
    
    
