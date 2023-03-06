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

        :param observed_length: length of the DNA sequence observed in the window.
        :type observed_length: int

        :param max_homopolymer_runs: maximum homopolymer runs.
        :type max_homopolymer_runs: int

        :param gc_range: range of GC content.
        :type gc_range: list

        :param undesired_motifs: undesired DNA motifs.
        :type undesired_motifs: list

        Example
            >>> from dsw import LocalBioFilter
            >>> bio_filter = LocalBioFilter(observed_length=8, \
                                            max_homopolymer_runs=2, gc_range=[0.4, 0.6], undesired_motifs=["GC"])
            >>> bio_filter.valid(dna_sequence="ACGTACGT")
            True
            >>> bio_filter.valid(dna_sequence="GCATGCAT")
            False
            >>> bio_filter.valid(dna_sequence="AAACCGGA")
            False

        .. notes::
            Reference [1] Nick Goldman et al. (2013) Nature

            Reference [2] Yaniv Erlich and Dina Zielinski (2017) Science

            Reference [3] William H. Press et al. (2020) Proceedings of the National Academy of Sciences

            Reference [4] Hannah F Lochel et al. (2021) Nucleic Acids Research

            If the maximum homopolymer runs (max_homopolymer_runs) is 1,
            "AA", "CC", "GG", "TT" cannot be included in tue valid DNA sequences.

            If the range of GC content (gc_range) is [0.4, 0.6],
            the GC content of valid DNA sequences must between 40% and 60%.

            If "GC" in the undesired DNA motifs (undesired_motifs), "GC" cannot be included in tue valid DNA sequences.
            This parameter could contain the restriction enzyme sites or some low compatibility DNA patterns.
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

        self.observed_length = observed_length
        self.max_homopolymer_runs = max_homopolymer_runs
        self.gc_range = gc_range
        self.undesired_motifs = undesired_motifs

    def valid(self, dna_sequence, only_last=True):
        """
        Judge whether the DNA sequence meets the local biochemical constraints.

        :param dna_sequence: DNA sequence to be judged.
        :type dna_sequence: str

        :param only_last: only check the DNA sequence of the last observed window.
        :type only_last: bool

        :return: judgement.
        :rtype: bool

        .. note::
            "only_last" parameter is used to save time.
            For most tree-based coding algorithms,
            it is not necessary to detect the sub DNA sequences observed in each window from scratch every time.
        """
        if only_last:
            observed_dna_sequence = dna_sequence[-self.observed_length:]
        else:
            observed_dna_sequence = dna_sequence

        for nucleotide in observed_dna_sequence:
            if nucleotide not in "ACGT":
                return False

        if self.max_homopolymer_runs is not None:
            for nucleotide in "ACGT":
                if nucleotide * (1 + self.max_homopolymer_runs) in observed_dna_sequence:
                    return False

        if self.undesired_motifs is not None:
            for special in self.undesired_motifs:
                if special in observed_dna_sequence:
                    return False
                reverse_complement = special.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
                reverse_complement = reverse_complement[::-1].upper()
                if reverse_complement in observed_dna_sequence:
                    return False

        if self.gc_range is not None:
            if len(observed_dna_sequence) >= self.observed_length:
                for index in range(len(observed_dna_sequence) - self.observed_length + 1):
                    sub_dna_sequence = observed_dna_sequence[index: index + self.observed_length]
                    gc_count = sub_dna_sequence.count("C") + sub_dna_sequence.count("G")
                    if gc_count > self.gc_range[1] * self.observed_length:
                        return False
                    if gc_count < self.gc_range[0] * self.observed_length:
                        return False
            else:
                gc_count = observed_dna_sequence.count("C") + observed_dna_sequence.count("G")
                if gc_count > self.gc_range[1] * self.observed_length:
                    return False
                at_count = observed_dna_sequence.count("A") + observed_dna_sequence.count("T")
                if at_count > (1 - self.gc_range[0]) * self.observed_length:
                    return False

        return True

    def __str__(self):
        info = self.screen_name + "\n"
        info += "maximum homopolymer runs : " + str(self.max_homopolymer_runs) + "\n"
        info += "local GC content range   : " + str(self.gc_range[0]) + " <= GC <= " + str(self.gc_range[1]) + "\n"
        info += "undesired DNA motifs     : " + str(self.undesired_motifs).replace("\"", "") + "\n"
        return info
