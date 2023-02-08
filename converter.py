import mosek
from exponential_rewriter import exponential_rewriter
from geometric_mean_rewriter import geometric_mean_rewriter

def main():
    instances = list()
    with open('benchmark-v2.txt', 'r') as file:
        for line in file:
            instances.append(line[0:-1]) # remove '\n'

    for instance in instances:
        exponential_rewriter("benchmark", instance, "benchmark_exp")
        geometric_mean_rewriter("benchmark", instance, "benchmark_gm")

    # instance = "ns1644855.mps.gz"
    # exponential_rewriter("benchmark", instance, "benchmark_exp")
    # geometric_mean_rewriter("benchmark", instance, "benchmark_gm")

main()