__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


import os
import math
import struct
import numpy
from dsw.algorithm import encode, decode
from dsw.monitor import Monitor


def read_bits_from_file(path, segment_length=120):
    monitor = Monitor()

    with open(path, mode="rb") as file:
        print("Read binary matrix from file: " + path)

        size = os.path.getsize(path)

        if segment_length > 0:
            # Set init storage matrix
            matrix = [[0 for _ in range(segment_length)] for _ in range(math.ceil(size * 8 / segment_length))]

            row = 0
            col = 0
            for byte_index in range(size):
                # Read a file as bytes
                one_byte = file.read(1)
                element = list(map(int, list(str(bin(struct.unpack("B", one_byte)[0]))[2:].zfill(8))))
                for bit_index in range(8):
                    matrix[row][col] = element[bit_index]
                    col += 1
                    if col == segment_length:
                        col = 0
                        row += 1
                monitor.output(byte_index + 1, size)
        else:
            matrix = []
            for byte_index in range(size):
                # Read a file as bytes
                one_byte = file.read(1)
                matrix += list(map(int, list(str(bin(struct.unpack("B", one_byte)[0]))[2:].zfill(8))))

            matrix = [matrix]

    return matrix, size * 8


def write_bits_to_file(path, matrix, bit_size):
    monitor = Monitor()

    with open(path, "wb+") as file:
        print("Write file from binary matrix: " + path)

        byte_size = int(bit_size / 8)
        # Change bit to byte (8 -> 1), and write a file as bytes
        bit_index = 0
        temp_byte = 0
        for row in range(len(matrix)):
            for col in range(len(matrix[0])):
                bit_index += 1
                temp_byte *= 2
                temp_byte += matrix[row][col]
                if bit_index == 8:
                    if byte_size > 0:
                        file.write(struct.pack("B", int(temp_byte)))
                        bit_index = 0
                        temp_byte = 0
                        byte_size -= 1
            monitor.output(row + 1, len(matrix))

    return True


def transcode(graph, start_index, matrix):
    monitor = Monitor()
    print("Encode by graph")
    oligos = []
    for index, bits in enumerate(matrix):
        oligos.append(encode(bits=bits, graph=graph, start_index=start_index))
        monitor.output(index + 1, len(matrix))

    print("Decode by graph")
    matrix = []
    for index, oligo in enumerate(oligos):
        matrix.append(decode(oligo=oligo, bit_length=120, graph=graph, start_index=start_index))
        monitor.output(index + 1, len(oligos))

    return matrix


if __name__ == "__main__":
    for screen_index in [1, 2, 3, 4, 5, 6]:
        vertices = numpy.load(file="../outputs/screen" + str(screen_index) + "[vertices].npy")
        si = numpy.random.choice(numpy.where(vertices == 1)[0])
        g = numpy.load(file="../outputs/screen" + str(screen_index) + "[graph].npy")

        m1, s = read_bits_from_file(path="../logo.svg", segment_length=120)
        m2 = transcode(graph=g, start_index=si, matrix=m1)
        write_bits_to_file(path="../outputs/logo" + str(screen_index) + ".svg", matrix=m2, bit_size=s)
        print()
