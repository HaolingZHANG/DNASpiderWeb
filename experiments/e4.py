from numpy import load, linspace, where
from os import path
from pickle import dump as psave
from experiments.code_repair import evaluate_repair_multiple_errors

task_seed = 2021
repeats = 10000
accessor, vertices = load("./results/data/a01[g].npy"), where(load("./results/data/a01[v].npy") == 1)[0]
dna_lengths = linspace(start=100, stop=800, num=8, dtype=int)
error_time = 4
for dna_length in dna_lengths:
    save_path = "./results/temp/multiple" + str(dna_length).zfill(4) + "." + str(error_time).zfill(2) + ".pkl"
    if not path.exists(save_path):
        records = evaluate_repair_multiple_errors(random_seed=task_seed, accessor=accessor, vertices=vertices,
                                                  observed_length=10, vt_length=10,
                                                  repeats=repeats, dna_length=dna_length,
                                                  error_times=error_time, check_iterations=error_time + 1)
        with open(save_path, "wb") as file:
            psave(records, file)
