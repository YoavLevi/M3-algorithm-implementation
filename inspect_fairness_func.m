function [alloc,load,max_load,stabilityness,balance,deviation] = inspect_fairness_func(q,M,lambda)
    min_alloc = floor(M.*q);
    alloc = min_alloc;
    rem = q - sum(min_alloc);
    for i = 1:rem
        add_mat = eye(numel(alloc));
        alloc_mat = alloc + add_mat;
        load = (alloc_mat.*lambda) ./ (q.*M);
        if sum(max(load') == min(max(load'))) > 1
            [~,chosen] = max(M.*(max(load') == min(max(load'))));
        else
            [~,chosen] = min(max(load'));
        end
        alloc(chosen) = alloc(chosen) + 1;
    end
    load = (alloc.*lambda) ./ (q.*M);
    max_load = max(load);
    stabilityness = min((M - alloc * lambda / q) ./ M);
    balance = std(load);
    if (rem && ~ismember(chosen,find(load==max(load))))
        error('unexpected error occured');
    end
    deviation = max(load - lambda);
end