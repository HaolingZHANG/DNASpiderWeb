__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


n_system = ["A", "C", "G", "T"]

b_system = [0, 1]


from dsw.biofilter import DefaultBioFilter, LocalBioFilter  # noqa: E402
from dsw.monitor import Monitor  # noqa: E402
from dsw.spectra import calculate_capacity  # noqa: E402
from dsw.spiderweb import find_vertices, connect_default_graph, connect_fixed_graph, connect_variable_graph  # noqa: E402
from dsw.coder import encode, decode  # noqa: E402
