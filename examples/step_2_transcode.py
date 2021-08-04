__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import os
import random
import numpy
from dsw.coder import encode, decode, obtain_fixed_length
from dsw.monitor import Monitor


def transcode(name, graph, start_indices, bit_length, test_time, fixed_length=None):
    """
    Run the encoding and decoding pipeline.

    :param name: name of this graph coder.
    :type name: str

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray

    :param start_indices: virtual vertices to start decoding (consistent with the encoding process).
    :type start_indices: list

    :param bit_length: length of bit array.
    :type bit_length: int

    :param test_time: transcoding test iteration.
    :type test_time: int

    :param fixed_length: length of oligo (it is required for fixed-length coder).
    :type fixed_length: int

    :raise: ValueError, if encode bit matrix is not equal to decoded bit matrix.

    :returns: decoded bit matrix and nucleotide number during encoding process.
    :rtype: tuple
    """
    random.seed(2021)
    encoded_matrix = [[random.randint(0, 1) for _ in range(bit_length)] for _ in range(test_time)]

    print("Evaluate " + name + ".")
    monitor, nucleotide_numbers = Monitor(), []
    for current, start_index in enumerate(start_indices):
        oligos, nucleotide_number = [], 0
        for bits in encoded_matrix:
            oligo = encode(bits=bits, graph=graph, start_index=start_index, length=fixed_length)
            nucleotide_number += len(oligo)
            oligos.append(oligo)

        nucleotide_numbers.append(nucleotide_number)

        decoded_matrix = []
        for oligo in oligos:
            decoded_matrix.append(decode(oligo=oligo, bit_length=bit_length, graph=graph, start_index=start_index))

        if encoded_matrix != decoded_matrix:
            raise ValueError("Encode matrix is not equal to decoded matrix!")

        monitor.output(current + 1, len(start_indices))

    return nucleotide_numbers


if __name__ == "__main__":
    l, t = 100, 100
    for index in [1, 2, 3, 4, 5, 6]:
        if not os.path.exists("../outputs/FLC-" + str(index) + " transcode.npy"):
            g = numpy.load(file="../entities/FLC" + str(index) + "[graph].npy")
            v = numpy.where(numpy.load(file="../entities/FLC" + str(index) + "[vertices].npy") == 1)[0]
            numbers = transcode(name="FLC-" + str(index), graph=g, start_indices=v, bit_length=l, test_time=t,
                                fixed_length=obtain_fixed_length(graph=g, name="FLC-" + str(index), bit_length=l))
            numbers = (10000 / numpy.array(numbers)) / 2
            numpy.save(file="../outputs/FLC-" + str(index) + " transcode.npy", arr=numbers)

        if not os.path.exists("../outputs/VLC-" + str(index) + " transcode.npy"):
            g = numpy.load(file="../entities/VLC" + str(index) + "[graph].npy")
            v = numpy.where(numpy.load(file="../entities/VLC" + str(index) + "[vertices].npy") == 1)[0]
            numbers = transcode(name="VLC-" + str(index), graph=g, start_indices=v, bit_length=l, test_time=t)
            numbers = (10000 / numpy.array(numbers)) / 2
            numpy.save(file="../outputs/VLC-" + str(index) + " transcode.npy", arr=numbers)
