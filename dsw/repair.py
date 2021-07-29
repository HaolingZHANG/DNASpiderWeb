__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import math
import numpy
from dsw import n_system


def repair_oligo(graph, oligo, start_index, need_logs):
    """
    Repair the oligo may contains one substitution error.

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray
    :param oligo: oligo waiting for recovery.
    :type oligo: str
    :param start_index: virtual vertex to start encoding.
    :type start_index: int
    :param need_logs: need to print log.
    :type need_logs: bool

    :return: repaired oligos (may contain multiple repair results).
    :rtype: list
    """
    vertex_index, index_queue = start_index, [start_index]
    for index, nucleotide in enumerate(oligo):
        try:
            used_indices = numpy.where(graph[vertex_index] >= 0)[0]
            [n_system[used_index] for used_index in used_indices].index(nucleotide)
            vertex_index = graph[vertex_index][n_system.index(nucleotide)]
        except ValueError:
            if need_logs:
                carbon_options, silicon_options = obtain_current(graph=graph, vertex_index=index_queue[-1])
                print("error occur in position " + str(index) + ", found " + oligo[index] +
                      " | " + str(carbon_options).replace("\'", "") + ", " + str(silicon_options))
            repaired_oligos = []

            # probe search
            if need_logs:
                carbon_options, silicon_options = obtain_current(graph=graph, vertex_index=index_queue[index])
                print("-" * 50)
                print("Probe (0) search error.")
                print("Check position = " + oligo[index] + " | " + str(index) + ", found " + oligo[index] +
                      " | " + str(carbon_options).replace("\'", "") + ", " + str(silicon_options))
            probe_oligos = check_repairability(graph=graph, oligo=oligo, index_queue=index_queue,
                                               occur_location=index, original_nucleotide=oligo[index],
                                               need_logs=need_logs)
            for probe_oligo in probe_oligos:
                if probe_oligo not in repaired_oligos:
                    repaired_oligos.append(probe_oligo)

            # recall search
            recall_queue = index_queue[::-1][:int(math.log(len(graph), len(n_system)))]
            for recall_index, original_vertex_index in enumerate(recall_queue):
                error = index - recall_index - 1
                if need_logs:
                    print("-" * 50)
                    print("Recall (" + str(-(recall_index + 1)) + ") search error.")
                    print("Check position = " + oligo[error] + " | " + str(error))
                recall_oligos = check_repairability(graph=graph, oligo=oligo, index_queue=index_queue,
                                                    occur_location=error, original_nucleotide=oligo[error],
                                                    need_logs=need_logs)
                for recall_oligo in recall_oligos:
                    if recall_oligo not in repaired_oligos:
                        repaired_oligos.append(recall_oligo)

            return repaired_oligos

        index_queue.append(vertex_index)

    return [oligo]


def check_repairability(graph, oligo, index_queue, occur_location, original_nucleotide, need_logs):
    """
    Perform saturation substitution at the selected position and obtain the oligos meeting repair conditions.

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray
    :param oligo: oligo waiting for saturation substitution in the specific location.
    :type oligo: str
    :param index_queue: queue of vertex indices (transcoding path in the graph).
    :type index_queue: list
    :param occur_location: the location that may occurring substitution (or specific location).
    :type occur_location: int
    :param original_nucleotide: original nucleotide in the mentioned-above location.
    :type original_nucleotide: str
    :param need_logs: need to print log.
    :type need_logs: bool

    :return: repaired oligos (may contain multiple repair results).
    :rtype: list
    """
    repair_oligos, observed_length = [], int(math.log(len(graph), len(n_system)))

    if occur_location + observed_length < len(oligo):
        check_segment = oligo[occur_location: occur_location + observed_length]
    else:
        check_segment = oligo[occur_location:]

    if need_logs:
        print("Original = [" + oligo + "]")
        if occur_location + observed_length < len(oligo):
            remain = len(oligo) - (occur_location + observed_length)
            print("Checking = [" + "*" * occur_location + check_segment + "*" * remain + "]")
        else:
            print("Checking = [" + "*" * occur_location + check_segment + "]")

    carbon_options, silicon_options = obtain_current(graph=graph, vertex_index=index_queue[occur_location])
    for replace_nucleotide in list(filter(lambda n: n != original_nucleotide, carbon_options)):
        segment = list(check_segment)
        segment[0] = replace_nucleotide
        segment = "".join(segment)
        if need_logs:
            print("Replace " + original_nucleotide + " to " + replace_nucleotide)
            if occur_location + observed_length < len(oligo):
                remain = len(oligo) - (occur_location + observed_length)
                print("         = [" + "*" * occur_location + segment + "*" * remain + "]" +
                      " | " + str(carbon_options).replace("\'", "") + ", " + str(silicon_options))
            else:
                print("         = [" + "*" * occur_location + segment + "]" +
                      " | " + str(carbon_options).replace("\'", "") + ", " + str(silicon_options))

        vertex_index, reliable = index_queue[occur_location], True
        for index, nucleotide in enumerate(segment):
            try:
                used_indices = numpy.where(graph[vertex_index] >= 0)[0]
                [n_system[used_index] for used_index in used_indices].index(nucleotide)
                if need_logs:
                    carbon_options, silicon_options = obtain_current(graph=graph, vertex_index=vertex_index)
                    if occur_location + observed_length < len(oligo):
                        remain = len(oligo) - (occur_location + observed_length)
                        print("           [" + "*" * (occur_location + index) + segment[index] +
                              "*" * (remain + len(segment) - index - 1) + "]" +
                              " | " + str(carbon_options).replace("\'", "") + ", " + str(silicon_options))
                    else:
                        print("           [" + "*" * (occur_location + index) + segment[index] + "]" +
                              " | " + str(carbon_options).replace("\'", "") + ", " + str(silicon_options))
                vertex_index = graph[vertex_index][n_system.index(nucleotide)]
            except ValueError:
                reliable = False

        if need_logs:
            print("Reliable? " + str(reliable))
        if reliable:
            obtained_oligo = list(oligo)
            obtained_oligo[occur_location] = replace_nucleotide
            obtained_oligo = "".join(obtained_oligo)
            repair_oligos.append(obtained_oligo)

    return repair_oligos


def obtain_current(graph, vertex_index):
    """
    Obtain carbon and silicon options based on current vertex index of the graph.

    :param graph: graph of DNA Spider-Web.
    :type graph: numpy.ndarray
    :param vertex_index: current vertex index.

    :returns: carbon and silicon options
    :rtype: tuple
    """
    used_indices = numpy.where(graph[vertex_index] >= 0)[0]
    carbon_options = [n_system[used_index] for used_index in used_indices]
    silicon_options = graph[vertex_index][~(graph[vertex_index] == -1)]
    return carbon_options, silicon_options
