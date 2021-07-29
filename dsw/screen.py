__author__ = "Zhang, Haoling [hlzchn@gmail.com]"


class DefaultScreen(object):

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
        raise NotImplementedError


class LocalScreen(DefaultScreen):

    def __init__(self, max_homopolymer_runs=None, max_gc_bias=None, ignore_motifs=None):
        """
        Initialize the screen of local biochemical constraints.

        :param max_homopolymer_runs: maximum homopolymer runs.
        If it is 1, "AA", "CC", "GG", "TT" cannot be included in tue valid oligos.
        :type max_homopolymer_runs: int
        :param max_gc_bias: maximum GC content bias.
        If it is 0.1, the GC content of valid oligos must between 40% and 60%.
        :type max_gc_bias: float
        :param ignore_motifs: DNA motifs to ignore, like restriction enzyme sites or low compatibility DNA segment.
        If "ACA" in this parameters,  "ACA" cannot be included in tue valid oligos.
        :type ignore_motifs: list
        """
        super().__init__(screen_name="Local")
        self.max_homopolymer_runs = max_homopolymer_runs
        self.max_gc_bias = max_gc_bias
        self.motifs = ignore_motifs

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

        if self.max_homopolymer_runs is not None:
            for nucleotide in "ACGT":
                if nucleotide * (1 + self.max_homopolymer_runs) in oligo:
                    return False

        if self.max_gc_bias is not None:
            if abs((float(oligo.count("C") + oligo.count("G")) / float(len(oligo))) - 0.5) > self.max_gc_bias:
                return False

        if self.motifs is not None:
            for special in self.motifs:
                if special in oligo:
                    return False
                reverse_complement = special.replace("A", "t").replace("C", "g").replace("G", "c").replace("T", "a")
                reverse_complement = reverse_complement[::-1].upper()
                if reverse_complement in oligo:
                    return False

        return True
