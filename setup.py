from setuptools import setup

setup(
    name="DNASpiderWeb",
    version="1.1",
    description="SPIDER-WEB generates coding algorithms with "
                "superior error tolerance and real-time information retrieval capacity",
    long_description="DNA has been considered a promising medium for storing digital information. "
                     "Previously, functions including bit-to-base transcoding and error correction are implemented by "
                     "independent algorithms. It would result in either increased computational complexity "
                     "or compromised error tolerance when attempting to correct various types of errors, "
                     "especially for insertions/deletions (indels). To address this issue, "
                     "we report a graph-based architecture, named SPIDER-WEB, providing an all-in-one coding solution "
                     "by generating customized algorithms with an efficient built-in error correction function. "
                     "SPIDER-WEB is able to correct a maximum of 4% edit errors in the DNA sequences "
                     "including substitution and indel, with only 5.5% logical redundancy. In addition, "
                     "SPIDER-WEB enables real-time retrieval of megabyte-level data with a 100x execution speed "
                     "improvement over conventional methods and hold the potential of practicability "
                     "for large-scale data storage applications at the exabyte-level.",
    author="Haoling ZHANG",
    author_email="zhanghaoling@genomics.cn",
    url="https://github.com/HaolingZHANG/DNASpiderWeb",
    packages=["dsw", "tests"],
    install_requires=["numpy", "networkx"],
    license="BGI-Research",
    classifiers=["Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
    keywords="DNA-based data storage, coding algorithm, automatic algorithm generator, probabilistic error correction",
)
