from setuptools import setup

setup(
    name="DNASpiderWeb",
    version="1.0",
    description="package for DNA Spider-Web",
    long_description="As a genetic material, DNA has become an attractive medium for storing digital information "
                     "gradually. In recent years, a growing number of functional biochemical operations and storage "
                     "environments were introduced, biochemical constraints are not confined to long single-nucleotide "
                     "repeats and abnormal GC content. However, trade-offs between information density and "
                     "compatibility to biochemical operations and algorithms require in-depth investigation, to create "
                     "transcoding algorithms considering novel biochemical constraints and their combinations. "
                     "Here, we design an automatic generator, named SPIDER-WEB, to create graph-based algorithms which "
                     "could be used directly or served as a benchmark for artificial algorithm design. "
                     "The generated coding algorithms provide high code rate, low parameter sensitivity, "
                     "certain self-repair ability, and valuable privacy protection potentiality. Through this work, "
                     "applicable algorithms with appropriate density-compatibility trade-off under arbitrary "
                     "local biochemical constraints could be generated in an automated way. "
                     "It is also suggested that more kinds of biochemical constraints can be further investigated "
                     "as more complex operations would be needed in future DNA storage systems.",
    author="Haoling ZHANG",
    author_email="zhanghaoling@genomics.cn",
    url="https://github.com/HaolingZHANG/DNASpiderWeb",
    packages=["dsw"],
    install_requires=["numpy", "scipy"],
    license="MIT",
    classifiers=["License :: OSI Approved :: MIT License",
                 "Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
    keywords="DNA Storage, Coding Algorithm, Automatic Algorithm Generator",
)
