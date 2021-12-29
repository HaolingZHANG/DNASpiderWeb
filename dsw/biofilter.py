__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


class DefaultBioFilter(object):

    def __init__(self, screen_name):
        """
        Initialize the default screen.

        :param screen_name: name of screen.
        :type screen_name: str
        """
        self.screen_name = screen_name

    def valid(self, dna_string):
        """
        Judge whether the DNA string meets the requirements.

        :param dna_string: DNA string to be judged.
        :type dna_string: str

        :raise: this interface needs to be implemented.

        :return: judgement.
        :rtype: bool
        """
        raise NotImplementedError("This interface \"def valid(dna_string)\" needs to be implemented.")


class LocalBioFilter(DefaultBioFilter):

    def __init__(self, max_homopolymer_runs=None, gc_range=None, undesired_motifs=None):
        """
        Initialize the screen of local biochemical constraints.

        :param max_homopolymer_runs: maximum homopolymer runs.
        :type max_homopolymer_runs: int

        :param gc_range: range of GC content.
        :type gc_range: list

        :param undesired_motifs: undesired DNA motifs.
        :type undesired_motifs: list

        ..notes:
            Reference [1] Nick Goldman et al. (2013) Nature

            Reference [2] Yaniv Erlich and Dina Zielinski (2017) Science

            Reference [3] William H. Press et al. (2020) Proceedings of the National Academy of Sciences

            Reference [4] Hannah F Löchel et al. (2021) Nucleic Acids Research

            If the maximum homopolymer runs (max_homopolymer_runs) is 1,
            "AA", "CC", "GG", "TT" cannot be included in tue valid DNA strings.

            If the range of GC content (gc_range) is [0.4, 0.6],
            the GC content of valid DNA strings must between 40% and 60%.

            If "GC" in the undesired DNA motifs (undesired_motifs), "GC" cannot be included in tue valid DNA strings.
            This parameter could contain the restriction enzyme sites or some low compatibility DNA patterns.

        ..example::
            >>> from dsw import LocalBioFilter
            >>> bio_filter = LocalBioFilter(max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["GC"])
            >>> bio_filter.valid(dna_string="ACGTACGT")
            True
            >>> bio_filter.valid(dna_string="GCATGCAT")
            False
            >>> bio_filter.valid(dna_string="AAACCGGA")
            False
        """
        super().__init__(screen_name="Local")
        self._max_homopolymer_runs = max_homopolymer_runs
        self._gc_range = gc_range
        self._undesired_motifs = undesired_motifs

    def valid(self, dna_string):
        """
        Judge whether the DNA string meets the local biochemical constraints.

        :param dna_string: DNA string to be judged.
        :type dna_string: str

        :return: judgement.
        :rtype: bool
        """
        for nucleotide in dna_string:
            if nucleotide not in "ACGT":
                return False

        if self._max_homopolymer_runs is not None:
            for nucleotide in "ACGT":
                if nucleotide * (1 + self._max_homopolymer_runs) in dna_string:
                    return False

        if self._gc_range is not None:
            gc = float(dna_string.count("C") + dna_string.count("G")) / float(len(dna_string))
            if gc < self._gc_range[0] or gc > self._gc_range[1]:
                return False

        if self._undesired_motifs is not None:
            for special in self._undesired_motifs:
                if special in dna_string:
                    return False
                reverse_complement = special.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
                reverse_complement = reverse_complement[::-1].upper()
                if reverse_complement in dna_string:
                    return False

        return True

    def __str__(self):
        info = self.screen_name + "\n"
        info += "maximum homopolymer runs : " + str(self._max_homopolymer_runs) + "\n"
        info += "local GC content range   : " + str(self._gc_range[0]) + "<= GC <=" + str(self._gc_range[1]) + "\n"
        info += "undesired DNA motifs     : " + str(self._undesired_motifs).replace("\"", "") + "\n"
        return info
