function num = return_num_req_keys(idx_on_ring,hash_values_sorted,hashed_keys)
    if (idx_on_ring==1)
        num = numel(find(hashed_keys > hash_values_sorted(end) & hashed_keys <= hash_values_sorted(idx_on_ring)));
    else
        num = sum(hashed_keys > hash_values_sorted(idx_on_ring-1) & hashed_keys <= hash_values_sorted(idx_on_ring));
    end
end

