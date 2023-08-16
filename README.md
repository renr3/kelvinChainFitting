# kelvinChainFitting

This repository is part of the conference paper entitled "An algorithm to compute the parameters of a Generalized Kelvin Chain model to represent aging creep of concrete", authored by Renan Rocha Ribeiro, Mara Monaliza Linhares Pereira, Valdirene Maria Silva Capuzzo, Rodrigo Lameiras, José Granja, Miguel Azenha, and presented on the conference XLIV Ibero-Latin American Congress on Computational Methods in Engineering (CILAMCE 2023), held in Porto from 13 to 16 November 2023.

The repository contains three modules related to viscoelastic modelling of concrete structures:

- creepEC2.py: a Python module that provides a class with methods to compute creep compliances accordingly to Eurocode 2 (CEN, “EN 1992-1-1 - Eurocode 2: Design of concrete structures - Part 1-1: General rules and rules for buildings,” 2004)
- creepMC2010.py: a Python module that provides a class with methods to compute creep compliances accordingly to fib Model Code 2010 (fib, fib Model Code for Concrete Structures 2010. Weinheim, Germany: Wiley-VCH Verlag GmbH & Co. KGaA, 2013. doi: 10.1002/9783433604090)
- fitKelvinChain.py: a Python module that provide a class with methods to compute the parameters of a Generalized Kelvin Chain model if a set of creep compliances is provided. Aging can also be handled. The Kelvin chain is derive by performing a non-linear curve fitting strategy, as suggested in the book "Mathematical modeling of creep and shrinkage of concrete" (1988), by  Z. P. Bazant.

A Jupyter notebook is included with some example codes on how to use the aforementioned modules. The codes generate the 2 figures provided in the conference paper. Also, two .xlsx files are included with data obtained from the software DIANA FEA and used to validate the Python modules for computing creep compliances.
