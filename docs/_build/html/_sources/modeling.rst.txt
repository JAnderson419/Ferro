Film Modeling Classes
======================

Introduction
-------------

This section details the modeling classes that exist in Ferro. Currently implemented is a hysteron-based multidomain Preisach model. This was implemented as a stepping stone towards a multidomain Landau model. LandauFilm is the base class for these modeled films, with the goal of eventually having a LandauSimple model for single-temperature measurements and a LandauFull model that handles the temperature-dependence of film hysteresis. LanduaDomain is a class that is used to represent individual ferroelectric domains in the film. 

I had hoped to implement the complete Landau model as part of my MS thesis work but did not get around to it. The functions are mostly implemented but I have not gotten around to testing them fully as I was having trouble getting good convergance on calculation of Landau Parameters from experimental data. If you are interested in working on this model to get it running, I would be glad to help you if you have any questions. 

LandauFilm
--------------
.. autoclass:: ferro.LandauFilm.LandauFilm
	:members:

LandauSimple
---------------
.. autoclass:: ferro.LandauFilm.LandauSimple
	:members:

LandauFull
------------
.. autoclass:: ferro.LandauFilm.LandauFull
	:members:

LandauDomain
------------
.. autoclass:: ferro.LandauFilm.LandauDomain
	:members:
