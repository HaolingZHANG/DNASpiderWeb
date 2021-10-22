__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


class DefaultBioFilter(object):

    def __init__(self, screen_name):
        """
        Initialize the default screen.

        :param screen_name: name of screen.
        :type screen_name:str
        """
        self.screen_name = screen_name

    def valid(self, oligo):
        """
        Judge whether the oligo meets the requirements.

        :param oligo: an oligo to be judged.
        :type oligo: str

        :raise: this interface needs to be implemented.

        :return: judgement.
        :rtype: bool
        """
        raise NotImplementedError("This interface \"def valid(oligo)\" needs to be implemented.")


class LocalBioFilter(DefaultBioFilter):

    def __init__(self, max_homopolymer_runs=None, gc_range=None, undesired_motifs=None):
        """
        Initialize the screen of local biochemical constraints.

        :param max_homopolymer_runs: maximum homopolymer runs.
        If it is 1, "AA", "CC", "GG", "TT" cannot be included in tue valid oligos.
        :type max_homopolymer_runs: int
        |cite| Nick Goldman et al. (2013) Nature

        :param gc_range: range of GC content.
        If it is [0.4, 0.6], the GC content of valid oligos must between 40% and 60%.
        :type gc_range: float
        |cite| Yaniv Erlich and Dina Zielinski (2017) Science

        :param undesired_motifs: undesired DNA motifs, like restriction enzyme sites or low compatibility DNA sequence.
        If "ACA" in this parameters,  "ACA" cannot be included in tue valid oligos.
        :type undesired_motifs: list
        """
        super().__init__(screen_name="Local")
        self._max_homopolymer_runs = max_homopolymer_runs
        self._gc_range = gc_range
        self._undesired_motifs = undesired_motifs

    def valid(self, oligo):
        """
        Judge whether the oligo meets the local biochemical constraints.

        :param oligo: an oligo to be judged.
        :type oligo: str

        :return: judgement.
        :rtype: bool
        """
        for nucleotide in oligo:
            if nucleotide not in "ACGT":
                return False

        if self._max_homopolymer_runs is not None:
            for nucleotide in "ACGT":
                if nucleotide * (1 + self._max_homopolymer_runs) in oligo:
                    return False

        if self._gc_range is not None:
            gc = float(oligo.count("C") + oligo.count("G")) / float(len(oligo))
            if gc < self._gc_range[0] or gc > self._gc_range[1]:
                return False

        if self._undesired_motifs is not None:
            for special in self._undesired_motifs:
                if special in oligo:
                    return False
                reverse_complement = special.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
                reverse_complement = reverse_complement[::-1].upper()
                if reverse_complement in oligo:
                    return False

        return True

    def __str__(self):
        info = self.screen_name + "\n"
        info += "maximum homopolymer runs : " + str(self._max_homopolymer_runs) + "\n"
        info += "local GC content range   : " + str(self._gc_range[0]) + "<= GC <=" + str(self._gc_range[1]) + "\n"
        info += "undesired DNA motifs     : " + str(self._undesired_motifs).replace("\"", "") + "\n"
        return info
