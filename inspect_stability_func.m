function [max_lambda_value] = inspect_stability_func(M,q)
    res = 1e-6;
    lambda = 0.5:res:1-res;
    orig_naive_func = sum(ceil(M.*q ./ lambda)-1) - q;
    results = orig_naive_func >= 0;
    indices = find(results);
    if isempty(indices)
        max_lambda_value = 0;
    else
        max_lambda_value = lambda(indices(end));
    end
end