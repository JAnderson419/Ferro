
Data Handling Classes
=================================

Introduction
-------------

This section includes documentation on the various data handling classes that are included in Ferro. The base class for all data classis is SampleData, which defines some common measurement parameters such as device dimensions and temperature. Derivced from this are LeakageData and HysteresisData. LeakageData is a container for I-V leakage measurement data while HysteresisData deals with P-V test data such as dynamic hysteresis, PUND, and FORC measurements.

Non-Class Functions
-----------------------
.. automodule:: ferro.HysteresisData
	:members:

SampleData
-----------
.. autoclass:: ferro.HysteresisData.SampleData
	:members:
	
HysteresisData
---------------
.. autoclass:: ferro.HysteresisData.HysteresisData
	:members:
	
LeakageData
--------------
.. autoclass:: ferro.HysteresisData.LeakageData
	:members: