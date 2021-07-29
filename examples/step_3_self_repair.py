__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import os
import random
import numpy
from dsw import b_system
from dsw.algorithm import encode
from dsw.repair import repair_oligo
from dsw.monitor import Monitor


def example():
    print("Case estimate the graph created by screening 1.")
    vertices = numpy.load(file="../outputs/screen1[vertices].npy")
    graph = numpy.load(file="../outputs/screen1[graph].npy")
    start_index = numpy.where(vertices == 1)[0][0]
    bits = b_system * 20
    wrong_location = 12
    right_oligo = encode(bits=bits, graph=graph, start_index=start_index)
    print("Right: " + right_oligo)
    print()
    for nucleotide in list(filter(lambda n: n != right_oligo[wrong_location], ["A", "C", "G", "T"])):
        wrong_oligo = list(right_oligo)
        wrong_oligo[wrong_location] = nucleotide
        wrong_oligo = "".join(wrong_oligo)
        print("Wrong: " + wrong_oligo)
        repaired_oligos = repair_oligo(graph=graph, oligo=wrong_oligo, start_index=start_index, need_logs=True)
        print("Repaired:")
        for repaired_oligo in repaired_oligos:
            print("       " + repaired_oligo + " | right? " + str(repaired_oligo == right_oligo))
        print()


def totally():
    random.seed(2021)
    bit_length, times, wrong_location = 30, 5, 12
    matrix = [[random.choice(b_system) for _ in range(bit_length)] for _ in range(times)]

    monitor = Monitor()
    for screen_index in [1, 2, 3, 4, 5, 6]:
        if not os.path.exists("../outputs/recover" + str(screen_index) + ".npy"):
            print("Totally estimate the graph created by screening " + str(screen_index) + ".")
            vertices = numpy.load(file="../outputs/screen" + str(screen_index) + "[vertices].npy")
            graph = numpy.load(file="../outputs/screen" + str(screen_index) + "[graph].npy")

            start_indices = numpy.where(vertices == 1)[0]
            results = numpy.zeros(shape=(len(start_indices), times * 3))
            for current, start_index in enumerate(start_indices):
                position = 0
                for array_index, bits in enumerate(matrix):
                    right_oligo = encode(bits=bits, graph=graph, start_index=start_index)
                    for nucleotide in list(filter(lambda n: n != right_oligo[wrong_location], ["A", "C", "G", "T"])):
                        wrong_oligo = list(right_oligo)
                        wrong_oligo[wrong_location] = nucleotide
                        wrong_oligo = "".join(wrong_oligo)
                        recover_oligos = repair_oligo(graph=graph, oligo=wrong_oligo, start_index=start_index,
                                                      need_logs=False)
                        found = right_oligo in recover_oligos
                        if found:
                            results[current][position] = len(recover_oligos)
                        else:
                            results[current][position] = 0
                        position += 1
                monitor.output(current + 1, len(start_indices))
            numpy.save(file="../outputs/recover" + str(screen_index) + ".npy", arr=results)


if __name__ == "__main__":
    example()
    totally()
