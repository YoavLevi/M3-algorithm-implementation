# M3-algorithm-implementation
A Matlab implementation of the M3 algorithm as described in [Beyond the Ring: Quantized Heterogeneous Consistent Hashing](https://icnp23.cs.ucr.edu/assets/papers/icnp23-final77.pdf). This repo also contains the code for the evaluations presented in the paper. 

# Documentation
## return_alloc(q,M)
An implementation function of the M3 algorithm.
- q - number of virtual servers in the system
- M - vector of service rates of the physical servers
The function returns an array indicating how many virtual servers should be assigned to each physical server according to the order presented in the service rates vector M.
## fnvhash_vec(msg)
A vectorized implementation function of the fnv-1a hash algorithm.
- msg - column vector of identifiers to hash
The function returns a column vector of hashed identifiers according to the order presented in the input vector.
## compare_algorithms_stability
A script file used to generate the output figure as presented in the evaluation section (Fig. 8).
# Usage
To run compare_algorithms_stability.m, you must first download the file ['unique_keys_hashed_154M.mat'](https://technionmail-my.sharepoint.com/:u:/g/personal/yoav1013_campus_technion_ac_il/ESr9spAiIn9GoH1CJtUmT2ABU9wn0YKuv8hdPD2byUHy6g?e=QmhfQR) into the same directory of the scripts.


