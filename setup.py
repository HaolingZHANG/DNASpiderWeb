from setuptools import setup

setup(
    name="DNASpiderWeb",
    version="1.0",
    description="package for Spider-Web",
    long_description="DNA has been considered a promising medium for storing digital information. "
                     "Despite the biochemical progress in DNA synthesis and sequencing, "
                     "novel coding algorithms need to be constructed under specific constraints. "
                     "As various biochemical constraints were introduced because of the "
                     "developing functional operations and storage carriers in recent years, "
                     "the existing methods may not fit additional constraints. "
                     "To solve this problem, we design a graph-based architecture, "
                     "named SPIDER-WEB, to generate coding algorithms under arbitrary local biochemical constraints. "
                     "These generated coding algorithms can serve state-of-the-art coding performances, "
                     "converting between digital data and DNA sequences or becoming a benchmark "
                     "for further artificial algorithm designs. "
                     "SPIDER-WEB also provides probabilistic error correction, "
                     "eliminating random errors from improper storage and reaction uncertainties "
                     "without multiple sequence alignment. Further, SPIDER-WEB can guarantee an ideal key strength "
                     "(equivalent to AES-256) under the existing DNA sequence design length.",
    author="Haoling ZHANG",
    author_email="zhanghaoling@genomics.cn",
    url="https://github.com/HaolingZHANG/DNASpiderWeb",
    packages=["dsw", "tests"],
    install_requires=["numpy", "scipy"],
    license="GPL",
    classifiers=["License :: OSI Approved :: GNU General Public License (GPL)",
                 "Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
    keywords="DNA-based data storage, coding algorithm, automatic algorithm generator, "
             "probabilistic error correction, encryption algorithm",
)
