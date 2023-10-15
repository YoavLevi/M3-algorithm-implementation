function hashes = fnvhash_vec( msg )

% FNV is for 32 bit int

FNV_PRIME_32 = 16777619;
FNV_OFFSET_32 = 2166136261;

n=size(msg,1);
hashes = uint32(zeros(n,1));
for j=1:n
    hash = FNV_OFFSET_32;
    for i = 1:length(msg(j,:))
       hash =  bitxor(hash,double(msg(j,i)));
%        hash =  mod((hash * FNV_PRIME_32),2^31-1);
       hash =  mod((hash * FNV_PRIME_32),2^32-1);
    end
   hashes(j) = hash;
end
end
