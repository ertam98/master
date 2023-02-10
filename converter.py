import mosek
from exponential_rewriter import exponential_rewriter
from geometric_mean_rewriter import geometric_mean_rewriter
from helpfunctions import streamprinter

def main():
    # instances = list()
    # with open('benchmark-v2.txt', 'r') as file:
    #     for line in file:
    #         instances.append(line[0:-1]) # remove '\n'

    # for instance in instances:
    #     exponential_rewriter("benchmark", instance, "benchmark_exp")
    #     geometric_mean_rewriter("benchmark", instance, "benchmark_gm")

    instance = "decomp2.mps.gz"
    exponential_rewriter("benchmark", instance, "debug_exp")
    geometric_mean_rewriter("benchmark", instance, "debug_gm")

    with mosek.Task() as task:
        task.set_Stream(mosek.streamtype.log, streamprinter)
        task.readdata("debug_exp/decomp2.ptf.gz")

        task.putdouparam(mosek.dparam.optimizer_max_time, 10.0)
        task.optimize()
        task.solutionsummary(mosek.streamtype.log)

main()