Modules
=======

.. image:: _static/logo.svg

**SPIDER-WEB** package consists of four modules:

- biochemical constraint module (`biofilter.py <https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/biofilter.py>`_): implementation of 'biochemical constraint filter'.

- transcoding process module (`spiderweb.py <https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/spiderweb.py>`_): implementation of 'coding algorithm generation', 'encoding process', 'decoding process', 'shuffled matrix creation for mapping shuffling', 'Varshamov-Tenengolts-based path check set', and 'obtained DNA string correction'.

- graph-based operation module (`graphized.py <https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/graphized.py>`_): implementation of 'data structure transformation between different graph representations', 'vertex search', 'path search', and 'capacity approximation'.

- fundamental operation module (`operation.py <https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/operation.py>`_): implementation of 'process monitor', 'large integer basic operation', and 'conversion between binary message, decimal number, and DNA string'.

.. toctree::
   :maxdepth: 2

Biochemical Constraint Module
------------------------------------------
.. autoclass:: dsw.biofilter.DefaultBioFilter
  :members:
  :undoc-members:
  :show-inheritance:
.. autoclass:: dsw.biofilter.LocalBioFilter
  :members:
  :undoc-members:
  :show-inheritance:

Transcoding Process Module
------------------------------------------
.. autofunction:: dsw.spiderweb.find_vertices
.. autofunction:: dsw.spiderweb.connect_valid_graph
.. autofunction:: dsw.spiderweb.connect_coding_graph
.. autofunction:: dsw.spiderweb.create_random_shuffles
.. autofunction:: dsw.spiderweb.encode
.. autofunction:: dsw.spiderweb.decode
.. autofunction:: dsw.spiderweb.set_vt
.. autofunction:: dsw.spiderweb.repair_dna
.. autofunction:: dsw.spiderweb.remove_nasty_arc

Graph-based Operation Module
------------------------------------------
.. autofunction:: dsw.graphized.approximate_capacity
.. autofunction:: dsw.graphized.path_matching
.. autofunction:: dsw.graphized.calculate_intersection_score
.. autofunction:: dsw.graphized.obtain_formers
.. autofunction:: dsw.graphized.obtain_latters
.. autofunction:: dsw.graphized.obtain_leaf_vertices
.. autofunction:: dsw.graphized.remove_useless
.. autofunction:: dsw.graphized.adjacency_matrix_to_accessor
.. autofunction:: dsw.graphized.accessor_to_adjacency_matrix
.. autofunction:: dsw.graphized.get_complete_accessor
.. autofunction:: dsw.graphized.latter_map_to_accessor
.. autofunction:: dsw.graphized.accessor_to_latter_map

Fundamental Operation Module
------------------------------------------
.. autoclass:: dsw.operation.Monitor
  :members:
  :undoc-members:
  :show-inheritance:
.. autofunction:: dsw.operation.calculus_addition
.. autofunction:: dsw.operation.calculus_subtraction
.. autofunction:: dsw.operation.calculus_multiplication
.. autofunction:: dsw.operation.calculus_division
.. autofunction:: dsw.operation.dna_to_number
.. autofunction:: dsw.operation.number_to_dna
.. autofunction:: dsw.operation.bit_to_number
.. autofunction:: dsw.operation.number_to_bit
