function [max_lambda_value] = inspect_real_stability_func_opt(M,real_load)
    lambda = 1;
    rho_array = (real_load.*lambda)./M;
    [max_rho] = max(rho_array);
    max_lambda_value = 1/max_rho;
end