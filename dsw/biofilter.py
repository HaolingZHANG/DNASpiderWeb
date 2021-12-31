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

    def __init__(self, observed_length, max_homopolymer_runs=None, gc_range=None, undesired_motifs=None):
        """
        Initialize the screen of local biochemical constraints.

        :param observed_length: length of the DNA string observed in the window.
        :type observed_length: int

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
            >>> bio_filter = LocalBioFilter(observed_length=8, \
                                            max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["GC"])
            >>> bio_filter.valid(dna_string="ACGTACGT")
            True
            >>> bio_filter.valid(dna_string="GCATGCAT")
            False
            >>> bio_filter.valid(dna_string="AAACCGGA")
            False
        """
        super().__init__(screen_name="Local")
        if max_homopolymer_runs is not None:
            if observed_length < max_homopolymer_runs:
                raise ValueError("The parameter \"observed_length\" must "
                                 + "longer than the parameter \"max_homopolymer_runs\"!")
        if undesired_motifs is not None:
            for index, undesired_motif in enumerate(undesired_motifs):
                if len(undesired_motif) > observed_length:
                    raise ValueError("The parameter \"observed_length\" must "
                                     + "longer than the length of any motif in the parameter \"undesired_motifs\"!")

        self._observed_length = observed_length
        self._max_homopolymer_runs = max_homopolymer_runs
        self._gc_range = gc_range
        self._undesired_motifs = undesired_motifs

    def valid(self, dna_string, only_last=True):
        """
        Judge whether the DNA string meets the local biochemical constraints.

        :param dna_string: DNA string to be judged.
        :type dna_string: str

        :param only_last: only check the DNA string of the last observed window.
        :type only_last: bool

        :return: judgement.
        :rtype: bool

        ..note::
            "only_last" parameter is interesting, which is used to save time.
            For most tree-based coding algorithms,
            it is not necessary to detect the sub DNA strings observed in each window from scratch every time.
        """
        if only_last:
            observed_dna_string = dna_string[-self._observed_length:]
        else:
            observed_dna_string = dna_string

        for nucleotide in observed_dna_string:
            if nucleotide not in "ACGT":
                return False

        if self._max_homopolymer_runs is not None:
            for nucleotide in "ACGT":
                if nucleotide * (1 + self._max_homopolymer_runs) in observed_dna_string:
                    return False

        if self._undesired_motifs is not None:
            for special in self._undesired_motifs:
                if special in observed_dna_string:
                    return False
                reverse_complement = special.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
                reverse_complement = reverse_complement[::-1].upper()
                if reverse_complement in observed_dna_string:
                    return False

        if self._gc_range is not None:
            if len(observed_dna_string) >= self._observed_length:
                for index in range(len(observed_dna_string) - self._observed_length + 1):
                    sub_dna_string = observed_dna_string[index: index + self._observed_length]
                    gc_count = sub_dna_string.count("C") + sub_dna_string.count("G")
                    if gc_count > self._gc_range[1] * self._observed_length:
                        return False
                    if gc_count < self._gc_range[0] * self._observed_length:
                        return False
            else:
                gc_count = observed_dna_string.count("C") + observed_dna_string.count("G")
                if gc_count > self._gc_range[1] * self._observed_length:
                    return False
                at_count = observed_dna_string.count("A") + observed_dna_string.count("T")
                if at_count > (1 - self._gc_range[0]) * self._observed_length:
                    return False

        return True

    def __str__(self):
        info = self.screen_name + "\n"
        info += "maximum homopolymer runs : " + str(self._max_homopolymer_runs) + "\n"
        info += "local GC content range   : " + str(self._gc_range[0]) + "<= GC <=" + str(self._gc_range[1]) + "\n"
        info += "undesired DNA motifs     : " + str(self._undesired_motifs).replace("\"", "") + "\n"
        return info
