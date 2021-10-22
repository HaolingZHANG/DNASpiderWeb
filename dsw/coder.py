__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from numpy import where
from dsw import n_system


def encode(bits, graph, start_index, length=None):
    """
    Encode a bit array by the specific graph of DNA Spider-Web.

    :param bits: bit array.
    :type bits: list

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray

    :param start_index: virtual vertex to start encoding.
    :type start_index: int

    :param length: length of oligo (it is required for fixed-length coder).
    :type length: int

    :return: an oligo encoded by this graph.
    :rtype: str
    """
    decimal_number, oligo, vertex_index = to_number(bits), "", start_index

    while decimal_number > 0:
        used_indices = where(graph[vertex_index] >= 0)[0]
        value = used_indices[int(decimal_number % len(used_indices))]
        nucleotide, vertex_index = n_system[value], graph[vertex_index][value]
        decimal_number //= len(used_indices)
        oligo += nucleotide

    if length is not None:
        while len(oligo) < length:
            value = where(graph[vertex_index] >= 0)[0][0]
            nucleotide, vertex_index = n_system[value], graph[vertex_index][value]
            oligo += nucleotide

    return oligo


def decode(oligo, bit_length, graph, start_index):
    """
    Decode an oligo by the specific graph of DNA Spider-Web.

    :param oligo: an oligo encoded bt this graph.
    :type oligo: str

    :param bit_length: length of the bit array.
    :type bit_length: int

    :param graph: graph of DNA Spider-Web (consistent with the encoding process).
    :type graph: numpy.ndarray

    :param start_index: virtual vertex to start decoding (consistent with the encoding process).
    :type start_index: int

    :return: a binary segment decoded by this graph.
    :rtype: list
    """
    decimal_number, array, vertex_index = 0, [], start_index
    for index, nucleotide in enumerate(oligo):
        used_indices = where(graph[vertex_index] >= 0)[0]
        value = [n_system[used_index] for used_index in used_indices].index(nucleotide)
        array.append((len(used_indices), value))
        vertex_index = graph[vertex_index][n_system.index(nucleotide)]

    for index, (out_degree, number) in enumerate(array[::-1]):
        decimal_number *= out_degree
        decimal_number += number

    bits = to_bits(decimal_number, bit_length)

    return bits


def to_number(bits):
    """
    Transform a bit array to the equivalent decimal number.

    :param bits: bit array.
    :type bits: list

    :return: equivalent decimal number of the inputted bit array.
    :rtype: int
    """
    value = 0
    for index, number in enumerate(bits[::-1]):
        value += (2 ** index) * number
    return int(value)


def to_bits(decimal_number, bit_length):
    """
    Transform a decimal number to the equivalent bit array with specific length.

    :param decimal_number: decimal number of the bit array.
    :param bit_length: default length of the bit array.

    :return: bit array.
    :rtype: list
    """
    def _to_list(number):
        values = []
        if number:
            values = _to_list(number // 2)
            return values + [number % 2]
        else:
            return values
    one_array = _to_list(decimal_number)
    return [0] * (bit_length - len(one_array)) + one_array
