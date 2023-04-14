# master
This is the code for my Master's thesis.

benchmark: The benchmark set from MIPLIB 2017 which contains 240 instances. (https://miplib.zib.de/)

benchmark-v2.txt: Text file with the namnes of all the instances in the benchmark set.

MyTask.py: mosek.Task() objected adapted for rewriting problems.

rewriter.py: Rewrites problems to geometric mean cone and exponential cone for analytical center (mps.gz/task.gz -> task.gz).

rewriter_presolve.py: Presolve problems (mps.gz/task.gz -> task.gz).

optimizer.py: Optimizes problems and exports stats to a csv-file.

calculate_stats.py: Calculates stats and create graphs.

helpfunctions.py: Miscellaneous functions.

findoptimalface.py: Solved problems to fix optimal face (mps.gz/task.gz -> task.gz).
