from setuptools import setup

setup(
    name="DNASpiderWeb",
    version="1.0",
    description="SPIDER-WEB generates coding algorithms with "
                "superior error tolerance and real-time information retrieval capacity",
    long_description="DNA has been considered a promising medium for storing digital information. "
                     "As an essential step in the DNA-based data storage workflow, "
                     "coding scheme is responsible to implement functions "
                     "including bit-to-base transcoding, error correction, etc. "
                     "In previous studies, these functions are normally realized by "
                     "introducing independent coding algorithms. "
                     "Here, we report a graph-based architecture, named SPIDER-WEB, "
                     "providing an all-in-one coding solution by generating customized algorithms automatically. "
                     "SPIDER-WEB supports correcting at maximum 4% edit errors in the DNA sequences "
                     "including substitution and insertion/deletion (indel), "
                     "while the required logical redundancy at only 5.5%. "
                     "Since no DNA sequence pretreatment is required for correcting and decoding processes, "
                     "SPIDER-WEB offers the function of real-time information retrieval, "
                     "which is 305.08 times faster than the speed of single-molecule sequencing techniques."
                     " Our retrieval process can improve 2 orders of magnitude faster compared to conventional one "
                     "under megabyte-level data and can be scalable to fit exabyte-level data. "
                     "Therefore, SPIDER-WEB holds the potential to improve the practicability "
                     "in large-scale data storage applications.",
    author="Haoling ZHANG",
    author_email="zhanghaoling@genomics.cn",
    url="https://github.com/HaolingZHANG/DNASpiderWeb",
    packages=["dsw", "tests"],
    install_requires=["numpy", "networkx"],
    license="GPL",
    classifiers=["License :: OSI Approved :: GNU General Public License (GPL)",
                 "Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
    keywords="DNA-based data storage, coding algorithm, automatic algorithm generator, probabilistic error correction",
)
