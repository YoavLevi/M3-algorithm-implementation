function [alloc] = return_alloc(q,M)
    min_alloc = floor(M.*q);
    alloc = min_alloc;
    rem = q - sum(min_alloc,2);
    for each_m = 1:size(M,1)
        for i = 1:rem(each_m)
            add_mat = eye(numel(alloc(each_m,:)));
            alloc_mat = alloc(each_m,:) + add_mat;
            load = (alloc_mat) ./ (q.*M(each_m,:));
            [~,chosen] = min(max(load'));
            alloc(each_m,chosen) = alloc(each_m,chosen) + 1;
        end
    end
end