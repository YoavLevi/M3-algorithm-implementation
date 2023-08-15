function [max_lambda_value] = inspect_real_stability_func_opt(M,real_load)
%     res = 1e-5;
%     lambda = 0:res:1-res;
    lambda = 1;
    rho_array = (real_load.*lambda)./M;
%     stability_array = M - real_load.*lambda;
%     results = all(stability_array > 0);
%     indices = find(results);
    [max_rho] = max(rho_array);
%     if isempty(indices)
%         max_lambda_value = 0;
%     else
%         max_lambda_value = lambda(indices(end));
%     end
    max_lambda_value = 1/max_rho;
end