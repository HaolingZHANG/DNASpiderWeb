__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import os
import math
import numpy
from dsw import n_system
from dsw.monitor import Monitor


def find_vertices(length, bio_filter, save_path=None):
    """
    Find default valid vertices based on the given constraints.

    :param length: length of the oligo in a vertex.
    :type length: int
    :param bio_filter: screening operation for identifying the valid oligo (required the given constraints).
    :type bio_filter: dsw.biofilter.DefaultBioFilter
    :param save_path: path to save file.
    :type save_path: str

    :return: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :rtype: numpy.ndarray
    """

    vertices, monitor = numpy.zeros(shape=(int(len(n_system) ** length),), dtype=numpy.int), Monitor()
    if save_path is not None and os.path.exists(save_path + "[vertices].npy"):
        vertices = numpy.load(file=save_path + "[vertices].npy")
    else:
        print("Find valid vertices in this length of oligo.")
        for vertex_index in range(len(vertices)):
            vertices[vertex_index] = bio_filter.valid(oligo=obtain_oligo(decimal_number=vertex_index, length=length))
            monitor.output(vertex_index + 1, len(vertices))

        if save_path is not None:
            numpy.save(file=save_path + "[vertices].npy", arr=vertices)

    valid_rate = numpy.sum(vertices) / len(vertices)

    if valid_rate > 0:
        print(str(round(valid_rate * 100, 2)) + "% (" + str(numpy.sum(vertices)) + ") valid vertices are found.")
        return vertices
    else:
        print("No valid vertex is found.")
        return None


def connect_default_graph(length, vertices, save_path=None):
    """
    Connect a default graph by valid vertices.

    :param length: length of the oligo in a vertex.
    :type length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param save_path: path to save file.
    :type save_path: str

    :return: default graph of DNA Spider-Web.
    The element of graph structure is:
    [latter index 1 (A), latter index 2 (C), latter index 2 (G), latter index 2 (T)]
    And the index of element is the current vertex index.
    :rtype: numpy.ndarray
    """
    if save_path is not None and os.path.exists(save_path + "[graph].npy"):
        transforms = numpy.load(file=save_path + "[graph].npy")
    else:
        valid_rate, monitor = numpy.sum(vertices) / len(vertices), Monitor()
        if valid_rate > 0:
            transforms = -numpy.ones(shape=(int(len(n_system) ** length), len(n_system)), dtype=numpy.int)
            for vertex_index in range(int(len(n_system) ** length)):
                if vertices[vertex_index]:
                    for position, latter_vertex_index in enumerate(obtain_latters(current=vertex_index, length=length)):
                        if vertices[latter_vertex_index]:
                            transforms[vertex_index][position] = latter_vertex_index

                monitor.output(vertex_index + 1, len(vertices))

            if save_path is not None:
                numpy.save(file=save_path + "[graph].npy", arr=transforms)
        else:
            print("No graph is created.")
            return None

    print("Default graph is created.")
    return transforms


def connect_fixed_graph(length, vertices, save_path=None):
    """
    Connect a fixed graph (fixed length code) by valid vertices.

    :param length: length of the oligo in a vertex.
    :type length: int

    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray

    :param save_path: path to save file.
    :type save_path: str

    :return: fixed graph of DNA Spider-Web.
    The element of graph structure is:
    [latter index 1 (A), latter index 2 (C), latter index 2 (G), latter index 2 (T)]
    And the index of element is the current vertex index.
    :rtype: numpy.ndarray
    """
    if save_path is not None and os.path.exists(save_path + "[graph].npy"):
        transforms = numpy.load(file=save_path + "[graph].npy")
    else:
        transforms, out_degree, monitor = None, None, Monitor()
        for out_degree in [4, 3, 2]:
            last_vertices = vertices.copy()
            print("Calculate " + str(out_degree) + " as minimum out-degree.")
            times = 1
            while True:
                print("Check the vertex collection requirement in round " + str(times) + ".")
                new_vertices = numpy.zeros(shape=(int(len(n_system) ** length),), dtype=numpy.int)
                saved_indices = numpy.where(last_vertices != 0)[0]
                for current, vertex_index in enumerate(saved_indices):
                    latter_indices = obtain_latters(current=vertex_index, length=length)
                    new_vertices[vertex_index] = numpy.sum(last_vertices[latter_indices]) >= out_degree
                    monitor.output(current + 1, len(saved_indices))

                changed = numpy.sum(last_vertices) - numpy.sum(new_vertices)
                rate = str(round(numpy.sum(new_vertices) / len(last_vertices) * 100, 2))
                print(rate + "% (" + str(numpy.sum(new_vertices)) + ") valid vertices are saved.")

                if not changed or numpy.sum(new_vertices) == 0:
                    break

                last_vertices = new_vertices
                times += 1

            last_vertices = new_vertices
            valid_rate, monitor = numpy.sum(last_vertices) / len(last_vertices), Monitor()
            if valid_rate > 0:
                vertices = last_vertices
                print("Preliminary connect graph.")
                transforms = -numpy.ones(shape=(int(len(n_system) ** length), len(n_system)), dtype=numpy.int)
                for vertex_index in range(int(len(n_system) ** length)):
                    if vertices[vertex_index]:
                        for position, latter_vertex_index in \
                                enumerate(obtain_latters(current=vertex_index, length=length)):
                            if vertices[latter_vertex_index]:
                                transforms[vertex_index][position] = latter_vertex_index

                    monitor.output(vertex_index + 1, len(vertices))
                break
            else:
                transforms = None

        if transforms is None:
            print("No graph is created.")
            return None

        print("Prune more arrows.")
        deleted_number = 0
        for vertex_index in range(len(transforms)):
            if len(numpy.where(transforms[vertex_index] >= 0)[0]) > out_degree:
                pruned_order = []
                for nucleotide in obtain_oligo(decimal_number=vertex_index, length=length)[::-1]:
                    position = n_system.index(nucleotide)
                    if position not in pruned_order:
                        pruned_order.append(position)
                if len(pruned_order) != 4:
                    for nucleotide in n_system:
                        position = n_system.index(nucleotide)
                        if position not in pruned_order:
                            pruned_order.append(position)
                for position in pruned_order:
                    transforms[vertex_index][position] = -1
                    deleted_number += 1
                    if len(numpy.where(transforms[vertex_index] >= 0)[0]) == out_degree:
                        break
            monitor.output(vertex_index + 1, len(vertices))
        print("Prune " + str(deleted_number) + " arrows.")

        times = 1
        while True:
            count = 0
            print("Check the vertex without former vertex in round " + str(times) + ".")
            saved_indices = numpy.where(vertices != 0)[0]
            for current, vertex_index in enumerate(saved_indices):
                former_indices = obtain_formers(current=vertex_index, length=length)
                if numpy.sum(vertices[former_indices]) == 0:
                    count += 1
                    vertices[vertex_index] = 0
                    transforms[vertex_index] = -1
                monitor.output(current + 1, len(saved_indices))

            print("Delete " + str(count) + " valid vertices without in-degree.")

            if not changed:
                break

            times += 1

        if save_path is not None:
            numpy.save(file=save_path + "[vertices].npy", arr=vertices)
            numpy.save(file=save_path + "[graph].npy", arr=transforms)

    print("Fixed graph is created.")
    return transforms


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

    :return: variable graph of DNA Spider-Web.
    The element of graph structure is:
    [latter index 1 (A), latter index 2 (C), latter index 2 (G), latter index 2 (T)]
    And the index of element is the current vertex index.
    :rtype: numpy.ndarray
    """
    if save_path is not None and os.path.exists(save_path + "[graph].npy"):
        transforms = numpy.load(file=save_path + "[graph].npy")
    else:
        times = 1
        while True:
            print("Check the vertex collection requirement in round " + str(times) + ".")
            new_vertices, monitor = numpy.zeros(shape=(int(len(n_system) ** length),), dtype=numpy.int), Monitor()
            saved_indices = numpy.where(vertices != 0)[0]
            for current, vertex_index in enumerate(saved_indices):
                latter_indices = obtain_latters(current=vertex_index, length=length)
                new_vertices[vertex_index] = numpy.sum(vertices[latter_indices]) > threshold
                monitor.output(current + 1, len(saved_indices))

            changed = numpy.sum(vertices) - numpy.sum(new_vertices)
            print(str(round(numpy.sum(new_vertices) / len(vertices) * 100, 2)) + "% (" + str(numpy.sum(new_vertices))
                  + ") valid vertices are saved.")

            if not changed or numpy.sum(new_vertices) == 0:
                break

            vertices = new_vertices
            times += 1

        valid_rate = numpy.sum(vertices) / len(vertices)
        if valid_rate > 0:
            transforms = -numpy.ones(shape=(int(len(n_system) ** length), len(n_system)), dtype=numpy.int)
            for vertex_index in range(int(len(n_system) ** length)):
                if vertices[vertex_index]:
                    for position, latter_vertex_index in enumerate(obtain_latters(current=vertex_index, length=length)):
                        if vertices[latter_vertex_index]:
                            transforms[vertex_index][position] = latter_vertex_index

                monitor.output(vertex_index + 1, len(vertices))

            if save_path is not None:
                numpy.save(file=save_path + "[vertices].npy", arr=vertices)
                numpy.save(file=save_path + "[graph].npy", arr=transforms)
        else:
            print("No graph is created.")
            return None

    print("Variable graph is created.")
    return transforms


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

    array = _to_list(decimal_number)
    return n_system[0] * (length - len(array)) + array


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
    if graph.shape[1] != len(n_system) or numpy.min(graph) < -1 or numpy.max(graph) > len(graph) - 1:
        raise ValueError("Wrong format in the graph of DNA Spider-Web")

    matrix, monitor = numpy.zeros(shape=(len(graph), len(graph)), dtype=numpy.int), Monitor()
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
    graph, monitor = -numpy.ones(shape=(len(matrix), len(n_system)), dtype=numpy.int), Monitor()
    for vertex_index, vertex in enumerate(matrix):
        next_indices = numpy.where(vertex == 1)[0].tolist()
        reference_latters = obtain_latters(current=vertex_index, length=int(math.log(len(matrix), len(n_system))))
        if list(set(next_indices) | set(reference_latters)) != reference_latters:
            raise ValueError("Wrong format in the adjacency matrix, "
                             + "which cannot be converted to equivalent compressed graph of DNA Spider-Web")
        saved_information = [index if index in next_indices else -1 for index in reference_latters]
        graph[vertex_index] = saved_information

        if need_log:
            monitor.output(vertex_index + 1, len(matrix))

    return graph
