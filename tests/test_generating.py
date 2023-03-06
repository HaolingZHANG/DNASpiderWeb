from numpy import array, all, where
from unittest import TestCase

from dsw import LocalBioFilter, find_vertices, connect_valid_graph, connect_coding_graph


class TestFindVertices(TestCase):

    def setUp(self):
        self.bio_filter = LocalBioFilter(observed_length=2, max_homopolymer_runs=2, gc_range=[0.5, 0.5])
        self.available_vertices = array([0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0], dtype=int)

    def test(self):
        vertices = find_vertices(observed_length=2, bio_filter=self.bio_filter).astype(int)
        self.assertEqual(all(self.available_vertices == vertices), True)


class TestGenerateValidGraph(TestCase):

    def setUp(self):
        self.vertices = array([0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0])
        self.valid_graph = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                                  [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                                  [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                                  [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])

    def test(self):
        graph = connect_valid_graph(observed_length=2, vertices=self.vertices)
        self.assertEqual(all(self.valid_graph == graph), True)


class TestGenerateCodingGraph(TestCase):

    def setUp(self):
        self.vertices = array([0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0])
        self.coding_graph = array([[-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1],
                                   [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                                   [-1, 1, 2, -1], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, 13, 14, -1],
                                   [-1, -1, -1, -1], [4, -1, -1, 7], [8, -1, -1, 11], [-1, -1, -1, -1]])

    def test(self):
        vertices, graph = connect_coding_graph(observed_length=2, vertices=self.vertices, threshold=1)
        self.assertEqual(all(where(self.vertices == 1)[0] == vertices.astype(int)), True)
        self.assertEqual(all(self.coding_graph == graph), True)
