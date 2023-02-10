# master
This is the code for my Master's thesis.

benchmark: The benchmark set from MIPLIB 2017 which contains 240 instances. (https://miplib.zib.de/)

benchmark-v2.txt: Text file with the namnes of all the instances in the benchmark set.

rewriter.py: Rewrites problems to geometric mean cone and exponential cone for analytical center (mps.gz -> ptf.gz).

optimizer.py: Optimizes problems and exports stats.

calculate_stats.py: Calculates stats and create graphs.
