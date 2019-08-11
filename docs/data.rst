
Data Handling Classes
=================================

Introduction
-------------

This section includes documentation on the various data handling classes that are included in Ferro. The base class for all data classis is SampleData, which defines some common measurement parameters such as device dimensions and temperature. Derivced from this are LeakageData and HysteresisData. LeakageData is a container for I-V leakage measurement data while HysteresisData deals with P-V test data such as dynamic hysteresis, PUND, and FORC measurements.

Non-Class Functions
-----------------------
.. automodule:: ferro.data
	:members:

SampleData
-----------
.. autoclass:: ferro.data.SampleData
	:members:
	
HysteresisData
---------------
.. autoclass:: ferro.data.HysteresisData
	:members:
	
LeakageData
--------------
.. autoclass:: ferro.data.LeakageData
	:members: