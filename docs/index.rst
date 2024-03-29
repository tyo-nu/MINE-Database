.. MINE-Database documentation master file, created by
   sphinx-quickstart on Wed Feb 10 14:27:07 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to MINE-Database's Documentation!
=========================================

Introduction
~~~~~~~~~~~~
MINE-Database, also referred to as Pickaxe, is a python library allows you to efficiently create reaction networks based on a set of reaction rules.

Some common use cases:

#. Predicting promiscuous enzymatic reactions in biological systems.
#. Searching for potential novel reaction pathways from starting compound(s) to target compound(s).
#. Annotating possible structures for unknown peaks in metabolomics datasets.
#. Predicting spontaneous chemical reactions which may be diverting flux from a pathway of interest.
#. Specifying custom reaction rules to extend reaction networks to include chemical reactions.

In all of these cases, you supply pickaxe with a set of starting compounds (as SMILES strings) and which set of reaction rules you would like to useand then Pickaxe does the rest.
Pickaxe creates a network expansion by applying these reaction rules iteratively to your starting set of compounds, going for as many generations as you specify.
There are many more advanced options and customizations you can add as well.

Getting Started
~~~~~~~~~~~~~~~
To get started, see :doc:`install`.

You can run pickaxe in two ways, in command-line mode (:doc:`command_line`) or using a template file (recommended) (:doc:`pickaxe_run`).
:doc:`pickaxe_run` also provides information about different compound filters you can apply to your pickaxe expansions.

For a list of inputs required for pickaxe, see :doc:`inputs`.

To learn how to create your own custom filters, see :doc:`custom_filters`.

An API reference is provided at :doc:`api_reference` if you need to see implementation details.

Finally, if you find yourself needing help or have feedback for us, please see :doc:`support`!

Contents
~~~~~~~~
.. toctree::
   :numbered:
   :maxdepth: 2

   install.rst
   command_line.rst
   pickaxe_run.rst
   inputs.rst
   custom_filters.rst
   thermodynamic_calculations.rst
   
   api_reference.rst

   support.rst