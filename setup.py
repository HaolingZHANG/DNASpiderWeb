from setuptools import setup

setup(
    name="DNASpiderWeb",
    version="1.0",
    description="package for Spider-Web",
    long_description="As a genetic material, DNA has become an attractive medium "
                     "for storing digital information gradually. "
                     "Besides the biochemical progress on DNA synthesis and sequencing, "
                     "novel coding algorithms need to be constructed catering to "
                     "the specific constraints in DNA storage. "
                     "In recent years, a growing number of functional biochemical operations "
                     "and storage environments were introduced, bringing in various biochemical constraints "
                     "including but not confined to long single-nucleotide repeats and abnormal GC content. "
                     "Given several local biochemical constraints and their combinations, "
                     "the code rate capacity and capacity-achieving coding algorithms require in-depth investigation. "
                     "In this work, we design an automatic generator, named SPIDER-WEB, "
                     "to create corresponding graph-based algorithms which could be used directly "
                     "or served as a benchmark for the construction of coding algorithms. "
                     "The main advantage of SPIDER-WEB is that "
                     "it provides an efficient way to applicable coding algorithms "
                     "for arbitrary local biochemical constraints in an automatic way "
                     "and support for probabilistic error correction and one-time pad encryption.",
    author="Haoling ZHANG",
    author_email="zhanghaoling@genomics.cn",
    url="https://github.com/HaolingZHANG/DNASpiderWeb",
    packages=["dsw", "tests"],
    install_requires=["numpy", "scipy"],
    license="GPL",
    classifiers=["License :: OSI Approved :: GNU General Public License (GPL)",
                 "Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
    keywords="DNA-based Storage, Coding Algorithm, Automatic Algorithm Generator, Probabilistic Error-Correction",
)
