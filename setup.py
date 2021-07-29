from setuptools import setup

setup(
    name="DNASpiderWeb",
    version="1.0",
    description="package for DNA Spider-Web",
    long_description="As a genetic material, DNA has become an attractive medium for storing digital information "
                     "gradually. Artificially-designed coding schemes may not be reliable under certain binary "
                     "patterns or biochemical constraints. Here, we report an automatic generator based on "
                     "Graph theory, named DNA Spider-Web, to establish end-to-end coding schemes with implicit rules "
                     "under the trade-off between information density and local feasibility (biochemical constraints: "
                     "homoploymer runs, GC content, and special motifs).",
    author="Haoling ZHANG",
    author_email="zhanghaoling@genomics.cn",
    url="https://github.com/HaolingZHANG/DNASpiderWeb",
    packages=["dsw"],
    install_requires=["numpy"],
    license="MIT",
    classifiers=["License :: OSI Approved :: MIT License",
                 "Programming Language :: Python :: 3",
                 "Operating System :: OS Independent"],
    keywords="DNA Storage, Coding Scheme",
)
