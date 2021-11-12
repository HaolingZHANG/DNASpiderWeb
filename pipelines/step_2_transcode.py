__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from os import path
from numpy import array, where, load, save
from random import seed, randint

from dsw import encode, decode
from dsw import Monitor


def transcode(graph, start_indices, task_seed, bit_length, test_time, fixed_length=None):
    seed(task_seed)
    encoded_matrix = [[randint(0, 1) for _ in range(bit_length)] for _ in range(test_time)]

    monitor, nucleotide_numbers = Monitor(), []
    for current, start_index in enumerate(start_indices):
        oligos, nucleotide_number = [], 0
        for bits in encoded_matrix:
            oligo = encode(bits=bits, graph=graph, start_index=start_index, length=fixed_length)
            nucleotide_number += len(oligo)
            oligos.append(oligo)

        nucleotide_numbers.append(nucleotide_number)

        # check the correctness of transcoding.
        decoded_matrix = []
        for oligo in oligos:
            decoded_matrix.append(decode(oligo=oligo, bit_length=bit_length, graph=graph, start_index=start_index))

        if encoded_matrix != decoded_matrix:
            raise ValueError("Encoded matrix is not equal to decoded matrix!")

        monitor.output(current + 1, len(start_indices))

    return nucleotide_numbers


if __name__ == "__main__":
    length, times = 100, 100
    for index in [1, 2, 3, 4, 5, 6, 7, 8]:
        if not path.exists("../results/VLC-" + str(index) + " transcode.npy"):
            print("Evaluate VLC-" + str(index) + ".")
            numbers = transcode(graph=load(file="../entities/VLC" + str(index) + "[g].npy"),
                                start_indices=where(load(file="../entities/VLC" + str(index) + "[v].npy") == 1)[0],
                                task_seed=2021, bit_length=length, test_time=times)
            numbers = ((length * times) / array(numbers)) / 2
            save(file="../results/VLC" + str(index) + "[t].npy", arr=numbers)
