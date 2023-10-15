function [max_lambda_value] = inspect_stability_real_ACH_func_opt(rl_ACH,normalized_mu)
    lambda = 1;
    valid_indices = find(max((rl_ACH'.*lambda)./(normalized_mu')) <= 1);
    rho_max = max((rl_ACH'.*lambda)./(normalized_mu'));
    max_lambda_value = 1/rho_max;
    
    
