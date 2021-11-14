__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


from os import path
from math import log, floor, ceil
from numpy import zeros, ones, array, linspace, sum, where, save, min, max

from dsw import n_system, Monitor


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
    decimal_number, values, vertex_index = 0, [], start_index
    for index, nucleotide in enumerate(oligo):
        used_indices = where(graph[vertex_index] >= 0)[0]
        value = [n_system[used_index] for used_index in used_indices].index(nucleotide)
        values.append((len(used_indices), value))
        vertex_index = graph[vertex_index][n_system.index(nucleotide)]

    for index, (out_degree, number) in enumerate(values[::-1]):
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


def find_vertices(length, bio_filter, save_path):
    """
    Find default valid vertices based on the given constraints.

    :param length: length of the oligo in a vertex.
    :type length: int
    :param bio_filter: screening operation for identifying the valid oligo (required the given constraints).
    :type bio_filter: dsw.biofilter.DefaultBioFilter
    :param save_path: path to save file.
    :type save_path: str
    """
    vertices, monitor = zeros(shape=(int(len(n_system) ** length),), dtype=int), Monitor()
    if not path.exists(save_path + "[v].npy"):
        print("Find valid vertices in this length of oligo.")
        for vertex_index in range(len(vertices)):
            vertices[vertex_index] = bio_filter.valid(oligo=obtain_oligo(decimal_number=vertex_index, length=length))
            monitor.output(vertex_index + 1, len(vertices))

        valid_rate = sum(vertices) / len(vertices)

        if valid_rate > 0:
            print(str(round(valid_rate * 100, 2)) + "% (" + str(sum(vertices)) + ") valid vertices are found.")
            save(file=save_path + "[v].npy", arr=vertices)


def connect_default_graph(length, vertices, save_path=None):
    """
    Connect a default graph by valid vertices.

    :param length: length of the oligo in a vertex.
    :type length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param save_path: path to save file.
    :type save_path: str
    """
    if not path.exists(save_path + "[g].npy"):
        print("Connect graph with valid vertices.")
        valid_rate, monitor = sum(vertices) / len(vertices), Monitor()
        if valid_rate > 0:
            transforms = -ones(shape=(int(len(n_system) ** length), len(n_system)), dtype=int)
            for vertex_index in range(int(len(n_system) ** length)):
                if vertices[vertex_index]:
                    for position, latter_vertex_index in enumerate(obtain_latters(current=vertex_index, length=length)):
                        if vertices[latter_vertex_index]:
                            transforms[vertex_index][position] = latter_vertex_index

                monitor.output(vertex_index + 1, len(vertices))

            if save_path is not None:
                save(file=save_path + "[g].npy", arr=transforms)
            print("Default graph is created.")
        else:
            print("No graph is created.")
    else:
        print("Default graph is created.")


def connect_fixed_graph(length, vertices, maximum_stride, save_path):
    """
    Connect a fixed graph (fixed length code) by valid vertices.

    :param length: length of the oligo in a vertex.
    :type length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param maximum_stride: maximum stride (pre transcoding step) in the graph.
    :type maximum_stride: int

    :param save_path: path to save file.
    :type save_path: str
    """
    assert length >= maximum_stride

    transforms, paths = None, []
    if not path.exists(save_path + "[g].npy"):
        monitor = Monitor()
        upper_bound, lower_bound, saved_information, saved_vertices = len(n_system), 2, None, None
        for stride in range(1, maximum_stride + 1):
            found_flag = False
            print("Calculate minimum out-degree with the stride " + str(stride) + " from "
                  + str(upper_bound) + " to " + str(lower_bound) + ".")
            for required_out_degree in list(range(max([2, lower_bound]),
                                                  min([4 ** stride + 1, upper_bound + 1])))[::-1]:
                print("Calculate " + str(required_out_degree) + " as minimum out-degree "
                      + "with the stride " + str(stride) + ".")
                last_vertices, times = vertices.copy(), 1
                while True:
                    print("Check the vertex collection requirement in round " + str(times) + ".")
                    saved_indices = where(last_vertices != 0)[0]
                    new_vertices = zeros(shape=(int(len(n_system) ** length),), dtype=int)
                    for current, vertex_index in enumerate(saved_indices):
                        level_indices = [vertex_index]
                        for _ in range(stride):
                            new_level_indices = []
                            for level_index in level_indices:
                                latter_indices = array(obtain_latters(current=level_index, length=length))
                                # check the midway vertices are also accessible.
                                access_indices = where(last_vertices[latter_indices] == 1)
                                new_level_indices += latter_indices[access_indices].tolist()
                            level_indices = new_level_indices
                        new_vertices[vertex_index] = sum(last_vertices[level_indices]) >= required_out_degree
                        monitor.output(current + 1, len(saved_indices))
                    last, current = sum(last_vertices), sum(new_vertices)
                    print("Last vertices have " + str(last) + ", and current " + str(current) + ".")
                    if sum(new_vertices) != 0:
                        if last != current:
                            last_vertices = new_vertices
                            times += 1
                        else:
                            found_flag = True
                            upper_bound = required_out_degree + 1
                            lower_bound = required_out_degree
                            saved_information = (required_out_degree, stride,
                                                 round(log(required_out_degree, 2) / stride, 3))
                            saved_vertices = last_vertices.copy()
                            paths.append([stride, required_out_degree])
                            print("Current information density = " + str(log(required_out_degree, 2) / stride))
                            print()
                            break
                    else:
                        print()
                        break
                if found_flag:
                    break
            if not found_flag:
                break

            upper_bound = upper_bound ** ((stride + 1) / stride)
            if upper_bound == floor(upper_bound):
                upper_bound = floor(upper_bound) - 1
            else:
                upper_bound = floor(upper_bound)

            lower_bound = lower_bound ** ((stride + 1) / stride)
            if lower_bound == ceil(lower_bound):
                lower_bound = ceil(lower_bound) + 1
            else:
                lower_bound = ceil(lower_bound)

            if upper_bound < lower_bound:
                break
            if saved_information[0] == len(n_system) ** stride:
                break

        if saved_information is not None and saved_vertices is not None:
            vertices = saved_vertices
            print("Generate graph with stride " + str(saved_information[1]) + ".")
            transforms = -ones(shape=(int(len(n_system) ** length), int(len(n_system) ** saved_information[1])),
                               dtype=int)
            used_vertices = where(vertices != 0)[0]
            for current, vertex_index in enumerate(used_vertices):
                level_indices, used_level_indices = [vertex_index], [vertex_index]
                for _ in range(saved_information[1]):
                    new_level_indices, new_used_level_indices = [], []
                    for level_index in level_indices:
                        latter_indices = array(obtain_latters(current=level_index, length=length))
                        new_level_indices += latter_indices.tolist()
                        access_indices = where(vertices[latter_indices] == 1)
                        new_used_level_indices += latter_indices[access_indices].tolist()
                    level_indices, used_level_indices = new_level_indices, new_used_level_indices
                level_indices, used_level_indices = array(level_indices), array(used_level_indices)

                selected_usages = (linspace(start=0, stop=len(used_level_indices) - 1,
                                            num=saved_information[0]) + 0.5).astype(int)
                used_level_indices = used_level_indices[selected_usages]
                locations = array([where(level_indices == used_level_index)[0][0]
                                   for used_level_index in used_level_indices])
                for location, leaf_vertex_index in zip(locations, used_level_indices):
                    transforms[vertex_index][location] = leaf_vertex_index
                monitor.output(current + 1, len(used_vertices))

            save(file=save_path + "[v].npy", arr=vertices)
            save(file=save_path + "[p].npy", arr=array(paths))
            save(file=save_path + "[g].npy", arr=transforms)
            print("Fixed graph is created.")
        else:
            print("No fixed graph is created.")
    else:
        print("Fixed graph is created.")


def connect_variable_graph(length, vertices, threshold, save_path=None):
    """
    Connect a variable graph (variable length code) by valid vertices and the threshold for minimum out-degree.

    :param length: length of the oligo in a vertex.
    :type length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param threshold: threshold for minimum out-degree.
    :type threshold: int

    :param save_path: path to save file.
    :type save_path: str
    """
    if not path.exists(save_path + "[g].npy"):
        times = 1
        while True:
            print("Check the vertex collection requirement in round " + str(times) + ".")
            new_vertices, monitor = zeros(shape=(int(len(n_system) ** length),), dtype=int), Monitor()
            saved_indices = where(vertices != 0)[0]
            for current, vertex_index in enumerate(saved_indices):
                latter_indices = obtain_latters(current=vertex_index, length=length)
                new_vertices[vertex_index] = sum(vertices[latter_indices]) > threshold
                monitor.output(current + 1, len(saved_indices))

            changed = sum(vertices) - sum(new_vertices)
            print(str(round(sum(new_vertices) / len(vertices) * 100, 2)) + "% (" + str(sum(new_vertices))
                  + ") valid vertices are saved.")

            if not changed or sum(new_vertices) == 0:
                break

            vertices = new_vertices
            times += 1

        valid_rate = sum(vertices) / len(vertices)
        if valid_rate > 0:
            transforms = -ones(shape=(int(len(n_system) ** length), len(n_system)), dtype=int)
            for vertex_index in range(int(len(n_system) ** length)):
                if vertices[vertex_index]:
                    for position, latter_vertex_index in enumerate(obtain_latters(current=vertex_index, length=length)):
                        if vertices[latter_vertex_index]:
                            transforms[vertex_index][position] = latter_vertex_index

                monitor.output(vertex_index + 1, len(vertices))

            save(file=save_path + "[v].npy", arr=vertices)
            save(file=save_path + "[g].npy", arr=transforms)
            print("Variable graph is created.")
        else:
            print("No variable graph is created.")
    else:
        print("Variable graph is created.")


def obtain_oligo(decimal_number, length):
    """
    Transform a decimal number to the equivalent oligo with specific length.

    :param decimal_number: decimal number of the oligo.
    :type decimal_number: int

    :param length: default length of the oligo.
    :type length: int

    :return: equivalent oligo of the decimal number.
    :rtype: str
    """
    def _to_list(number):
        values = ""
        if number:
            values = _to_list(number // len(n_system))
            return values + n_system[number % len(n_system)]
        else:
            return values

    one_array = _to_list(decimal_number)
    return n_system[0] * (length - len(one_array)) + one_array


def obtain_number(oligo):
    """
    Transform an oligo to the equivalent decimal number.

    :param oligo: required oligo.
    :type oligo: str

    :return: equivalent decimal number of the inputted oligo.
    :rtype: int
    """
    values, value = list(map(n_system.index, oligo)), 0
    for index, number in enumerate(values[::-1]):
        value += (len(n_system) ** index) * number
    return int(value)


def obtain_formers(current, length):
    """
    Obtain former vertex indices based on the current vertex index.

    :param current: current vertex index.
    :type current: int

    :param length: length of the oligo in a vertex.
    :type length: int

    :return: former vertex indices.
    :rtype: list
    """
    formers = []
    for former_value in range(len(n_system)):
        former = current // len(n_system) + former_value * int(len(n_system) ** (length - 1))
        formers.append(former)
    return formers


def obtain_latters(current, length):
    """
    Obtain latter vertex indices based on the current vertex index.

    :param current: current vertex index.
    :type current: int

    :param length: length of the oligo in a vertex.
    :type length: int

    :return: latter vertex indices.
    :rtype: list
    """
    latters = []
    for latter_value in range(len(n_system)):
        latter = int((current * len(n_system) + latter_value) % (len(n_system) ** length))
        latters.append(latter)
    return latters


def to_adjacency_matrix(graph, maximum_length=8, need_log=False):
    """
    Transform the graph of DNA Spider-Web (compressed) to an adjacency matrix of the uncompressed graph.

    :param graph: graph of DNA Spider-Web.
    :type: numpy.ndarray

    :param maximum_length: maximum vertex length (like 8 in general) of the adjacency matrix (size is (4^l)^2).
    :type maximum_length: int

    :param need_log: need to print log.
    :type need_log: bool

    :raise MemoryError: when you generate a large adjacency matrix that your memory cannot allocate.
    :raise ValueError: when you input a graph of DNA Spider-Web with wrong format.

    :return: adjacency matrix of the uncompressed graph.
    :rtype: numpy.ndarray
    """
    if len(graph) >= 4 ** maximum_length:
        raise MemoryError("Unable to allocate adjacency matrix when length of oligo (vertex) is more than 7.")
    if graph.shape[1] != len(n_system) or min(graph) < -1 or max(graph) > len(graph) - 1:
        raise ValueError("Wrong format in the graph of DNA Spider-Web")

    matrix, monitor = zeros(shape=(len(graph), len(graph)), dtype=int), Monitor()
    for vertex_index, vertex in enumerate(graph):
        matrix[vertex_index][vertex[vertex >= 0]] = 1

        if need_log:
            monitor.output(vertex_index + 1, len(graph))

    return matrix


def to_graph(matrix, need_log=False):
    """
    Transform the adjacency matrix of the uncompressed graph to a graph of DNA Spider-Web (compressed).

    :param matrix: adjacency matrix of the uncompressed graph of DNA Spider-Web.
    :type matrix: numpy.ndarray

    :param need_log: need to print log.
    :type need_log: bool

    :raise ValueError: when you input an adjacency matrix with wrong format.

    :return: equivalent (compressed) graph of DNA Spider-Web. (compress rate = len(n_system) / len(matrix))
    :rtype: numpy.ndarray
    """
    graph, monitor = -ones(shape=(len(matrix), len(n_system)), dtype=int), Monitor()
    for vertex_index, vertex in enumerate(matrix):
        next_indices = where(vertex == 1)[0].tolist()
        reference_latters = obtain_latters(current=vertex_index, length=int(log(len(matrix), len(n_system))))
        if list(set(next_indices) | set(reference_latters)) != reference_latters:
            raise ValueError("Wrong format in the adjacency matrix, "
                             + "which cannot be converted to equivalent compressed graph of DNA Spider-Web")
        saved_information = [index if index in next_indices else -1 for index in reference_latters]
        graph[vertex_index] = saved_information

        if need_log:
            monitor.output(vertex_index + 1, len(matrix))

    return graph
