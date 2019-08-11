AixACCT Data Import
======================

Introduction
-------------
Ferro implements functions to directly convert
data exported from an aixACCT brand ferroelectric
tester as a .dat file into the SampleData object
format used by the package. This is currently implemented
for hysteresis and leakage current measurements, although
it has been implemented in a way that other measurement types
could easily be imported if SampleData instances were created
for them.



In order to read aixACCT TF-x000 tester data into Ferro, it must
first be exported as text using the "Export as ASCII" file menu option.
This will save a .dat file to disk, which can be parsed by the ferro commands
in the aixACCT module as follows::


    from ferro import aixacct

    my_tf_data = r'/Path/to/my/file/here.dat'
    data_dict = aixacct.read_tfdata(my_tf_data)
    loaded_SampleData_obj_list = aixacct.load_tfdata(data_dict)

The data_dict in the above code contains all of the information contained
in the original .dat file, excluding the summary table at the top of the document.
The structure of the dictionary is as follows:

::

    'here.dat':
        'meastype': MeasEnum.HYSTERESIS,
        'datatables':
            'Table 1':
                'metadata':
                    Waveform: triangle
                    SampleName: RTWhiteB
                    Area [mm2]: 0.01
                    Thickness [nm]: 255
                    ...
                'data':
                    '0.000000e+000\t3.235215e-004\t ...'
            'Table 2':
                'metadata':
                    Waveform: triangle
                    SampleName: RTWhiteB
                    Area [mm2]: 0.01
                    Thickness [nm]: 255
                    ...
                'data':
                    '0.000000e+000\t3.235215e-004\t ...'
                ...

One can query data in this structure using standard python dictionary syntax such as::

    >>>data_dict['here.dat']['datatables']['Table 1']['metadata']['SampleName']
    'RTWhiteB'

The load_tfdata function then takes this dict and loads the relevant parameters from
each present data table into a corresponding SampleData object and returns all tables
from the file as a list. From here, one can access individual data tables in the list
for analysis or loop over the list to perform operations on each measurement.


MeasEnum
-----------
.. autoclass:: ferro.aixacct.MeasEnum
	:members:


Non-Class Functions
-----------------------
.. automodule:: ferro.aixacct
	:members: