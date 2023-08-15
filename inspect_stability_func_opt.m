function [max_lambda_value] = inspect_stability_func_opt(M,q_i)
%     res = 1e-6;
%     lambda = 0:res:1-res;
    lambda = 1;
    rho_array = (q_i .* lambda) ./ (M.* sum(q_i));
    [max_rho] = max(rho_array);
    max_lambda_value = 1/max_rho;
end