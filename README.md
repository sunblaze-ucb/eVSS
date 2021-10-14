# submission-eVSS
## Note to AE reviewers
The trusted setup version will be slower than reported numbers since it run N verifications, so the actual running time of our binary is Prover's time + N * (verifier's time).

To verify our time measurement on trusted setup version, open https://github.com/sunblaze-ucb/eVSS/blob/main/KZG_ext/src/kzg_main.cpp and check line 29-31

To verify our time measurement on Virgo version, open https://github.com/sunblaze-ucb/eVSS/blob/main/Virgo/src/linear_gkr/zk_prover.cpp and make sure all functions with name "sumcheck" as a substring has been covered by std::chrono::high_resolution_clock

