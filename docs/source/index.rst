Documentation of SPIDER-WEB
===========================

.. image:: _static/logo.svg

**SPIDER-WEB** is a graph-based architecture,
which provides tools for creating stable, repairable, and encryptible algorithms
under arbitrary local biochemical constraints.

First of all, as a generator, **SPIDER-WEB** provides graph-based coding algorithms
under arbitrary local biochemical constraint combinations.

.. image:: _static/strategy_1.png

Afterwards, as a corrector,
when one or more edit (including insertion, deletion, or substitution) errors occur during the DNA sequencing process,
**SPIDER-WEB** supplies a probabilistic repair strategy.

.. image:: _static/strategy_2.png

Finally, as a confounder, **SPIDER-WEB** offers certain encryption capabilities to against eavesdropping.

.. image:: _static/strategy_3.png

.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   structure
   modules
   customization
   reference
