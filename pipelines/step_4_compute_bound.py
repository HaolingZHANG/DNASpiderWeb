__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import numpy
from dsw.evaluator import haemers


if __name__ == "__main__":
    for screen_index in [1, 2, 3, 4, 5, 6]:
        graph = numpy.load(file="../outputs/VLC" + str(screen_index) + "[graph].npy")
        theta = haemers(graph=graph)
