Customization
=============

.. image:: _static/logo.svg

You can create your customized biochemical constraint filter (as the biochemical constraints set)
by inheriting `DefaultBioFilter <https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/biofilter.py#L4>`_.
For example:

.. code-block:: python

    from dsw import DefaultBioFilter

    class RegionalizedGCFilter(DefaultBioFilter):

        def __init__(self, window_length, gc_bias):
            super().__init__(screen_name="Regionalized GC content constraint")
            self._window_length = window_length
            self._gc_bias = gc_bias  # bias based on 0.5 (or 50%)

        def valid(self, dna_string):
            if len(dna_string) >= self._window_length:
                for index in range(len(dna_string) - self._window_length + 1):  # judge in a window
                    regional_dna_string = dna_string[index: index + self._window_length]
                    gc_count = regional_dna_string.count("C") + regional_dna_string.count("G")
                    if gc_count > (0.5 + self._gc_bias) * self._window_length:
                        return False
                    if gc_count < (0.5 - self._gc_bias) * self._window_length:
                        return False

            else:
                gc_count = dna_string.count("C") + dna_string.count("G")
                if gc_count > (0.5 + self._gc_bias) * self._window_length:
                    return False
                at_count = dna_string.count("A") + dna_string.count("T")
                if at_count > (0.5 + self._gc_bias) * self._window_length:
                    return False

            return True

Here is an investigated example in this work named `LocalBioFilter <https://github.com/HaolingZHANG/DNASpiderWeb/blob/main/dsw/biofilter.py#L30>`_.

Besides, according to **SPIDER-WEB**, you can obtain the corresponding variable-length algorithms
based on the customized biochemical constraints, for example:

.. code-block:: python

    from dsw import LocalBioFilter, find_vertices, connect_coding_graph

    bio_filter = LocalBioFilter(max_homopolymer_runs=4, gc_range=[0.1, 0.3], undesired_motifs=["GCC"])
    vertices = find_vertices(observed_length=10, bio_filter=bio_filter)
    accessor = connect_coding_graph(observed_length=10, vertices=vertices, threshold=1)
