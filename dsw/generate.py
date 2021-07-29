__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import os
import numpy
from dsw import n_system
from dsw.monitor import Monitor


def find_vertices(length, screen, save_path=None):
    """
    Find original valid vertices based on the given constraints.

    :param length: length of the oligo in a vertex.
    :type length: int
    :param screen: screening operation for identifying the valid oligo (required the given constraints).
    :type screen: dsw.screen.DefaultScreen
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
            vertices[vertex_index] = screen.valid(oligo=obtain_oligo(decimal_number=vertex_index, length=length))
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


def connect_graph(length, vertices, threshold, save_path=None):
    """
    Connect graph by valid vertices and the threshold for minimum out-degree.

    :param length: length of the oligo in a vertex.
    :type length: int
    :param vertices: vertex accessor, in each cell, True is valid vertex and False is invalid vertex.
    :type vertices: numpy.ndarray
    :param threshold: threshold for minimum out-degree.
    :type threshold: int
    :param save_path: path to save file.
    :type save_path: str

    :return: graph of DNA Spider-Web.
    :rtype: numpy.ndarray
    """
    if save_path is not None and os.path.exists(save_path + "[graph].npy"):
        vertices = numpy.load(file=save_path + "[vertices].npy")
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
            print(str(round(numpy.sum(new_vertices) / len(vertices) * 100, 2)) + "% (" + str(numpy.sum(new_vertices)) +
                  ") valid vertices are saved.")

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

    print("Graph with " + str(numpy.sum(vertices)) + " vertices is created.")
    return transforms


def obtain_oligo(decimal_number, length):
    """
    Transform a decimal number to the equivalent oligo with specific length.

    :param decimal_number: decimal number of the oligo.
    :param length:  original length of the oligo.

    :return: equivalent oligo of the decimal number.
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
