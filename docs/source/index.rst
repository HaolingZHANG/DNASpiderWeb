Documentation of SPIDER-WEB
===========================

.. image:: _static/logo.svg

**SPIDER-WEB** is a graph-based architecture in DNA-based storage,
which provides tools for creating stable, repairable, and encryptible algorithms
under arbitrary local biochemical constraints.

First of all, as an algorithm generator, **SPIDER-WEB** provides **graph-based encoding**
under arbitrary local biochemical constraint combinations.

.. image:: _static/strategy_1.svg

Afterwards, as a decoding repairer, **SPIDER-WEB** supplies **path-based correcting**
for obtained DNA strings with one or more edit errors (including substitution, insertion, and deletion) probabilistically.

.. image:: _static/strategy_2.svg

Last but not least, as a data confounder, ****SPIDER-WEB**** offers **mapping shuffling**,
to strengthen the capability of privacy protection.

.. image:: _static/strategy_3.svg


.. toctree::
   :maxdepth: 1
   :hidden:

   installation
   structure
   modules
   customization
   reference
