
import re
import numpy as np
import data as hd
from os.path import basename
from enum import Enum

class MeasEnum(Enum):
    HYSTERESIS = 1
    FATIGUE = 2
    PULSE = 3
    LEAKAGE = 4
meas_struct = {
    MeasEnum.HYSTERESIS: {
        'datatype': hd.HysteresisData,
        'metadata': {
            'area': 'Area [mm2]',
            'freq': 'Hysteresis Frequency [Hz]',
            'thickness': 'Thickness [nm]'
        },
        'datatable': {
            'time': 0,
            'voltage': 1,
            'current': 3,
            'polarization': 4
        },
        'multiplier': {
            'thickness': 1E-7,
            'polarization': 1E-6,
            'area': 1E-2
        }
    },
    MeasEnum.FATIGUE: None,
    MeasEnum.PULSE: None,
    MeasEnum.LEAKAGE: {
        'datatype': hd.LeakageData,
        'metadata': {
            'area': 'Area [mm2]',
            'thickness': 'Thickness [nm]'
        },
        'datatable': {
            'lcm_voltage': 0,
            'lcm_current': 1,
        },
        'multiplier': {
            'thickness': 1E-7,
            'area': 1E-2,
            'lcm_current': 1E-6
        }
    }

}

def check_datatype(filepath):
    f = open(filepath, encoding='cp1252')
    firstline = f.readline()
    f.close()
    if re.match('DynamicHysteresisResult', firstline):
        return MeasEnum.HYSTERESIS
    elif re.match('Fatigue', firstline):
        return MeasEnum.FATIGUE
    elif re.match('PulseResult', firstline):
        return MeasEnum.PULSE
    elif re.match('LeakageResult', firstline):
        return MeasEnum.LEAKAGE

def check_istable(line, datatype):
    '''

    Parameters
    ----------
    line
    datatype

    Returns
    -------

    '''
    if datatype is MeasEnum.LEAKAGE:
        m = re.match(r'^Voltage \[V\]', line)
    else:
        m = re.match(r'^Time \[s\]', line)
    return bool(m)

def read_tfdata(filepath):
    '''
    Read an AixACCT .dat text file

    Returns data in a dictionary of following structure :
        filename:
            'table #':
                'metadata':
                    Waveform: triangle
                    SampleName: RTWhiteB
                    Area [mm2]: 0.01
                    Thickness [nm]: 255
                    ...
               'data':
                    '0.000000e+000\t3.235215e-004\t ...'


    Parameters
    ----------
    filepath: str
        Path (inc. filename) of data to parse

    Returns
    -------
    table_dict: dict
        Dictionary containing parsed text file.
    '''
    summary_table_read = False
    is_data_table = False
    is_data_table_header = False
    is_global_header = True
    datatype = check_datatype(filepath)
    global_metadata = {}
    table_metadata = {}
    datastr = []
    tableheaderstr = ''
    table_key = ''

    filekey = basename(filepath).split('.')[0]
    table_dict= {
        filekey: {
            'meastype': datatype,
            'datatables': {}
        }
    }
    with open(filepath, encoding='cp1252') as f:
        for line in f:
            line.rstrip(r'\n')
            if re.match('^(DynamicHysteresis|Data Table|Pulse|Leakage)$', line):
                summary_table_read = True
            if not summary_table_read:
                continue
            else:
                if is_global_header:
                    if re.match(r'Table (\[*\d*.*\d*\]*)', line):
                        is_global_header = False
                        is_data_table_header = True
                        table_key = line
                    m = re.match(r'(.*): (.*)', line)
                    if m:
                        global_metadata[m.group(1)] = m.group(2)
                elif is_data_table_header:
                    if check_istable(line, datatype):
                        is_data_table = True
                        is_data_table_header = False
                        tableheaderstr = line
                        continue
                    else:
                        m = re.match(r'(.*): (.*)', line)
                        if m:
                            table_metadata[m.group(1)] = m.group(2)
                elif is_data_table:
                    if re.match(r'^(\n)?$', line):
                        table_dict[filekey]['datatables'].update({
                            table_key: {
                                'metadata': {**global_metadata,
                                             **table_metadata},
                                'dataheader': tableheaderstr.split('\t'),
                                'data': datastr
                            }
                        }
                        )
                    elif re.match(r'Table (\[*\d*.*\d*\]*)', line):
                        table_key = line
                        is_data_table_header = True
                        is_data_table = False
                        table_metadata.clear()
                        datastr = []
                    else:
                        datastr.append(line)
    # add last table in file to dict
    table_dict[filekey]['datatables'].update({
        table_key: {
            'metadata': {**global_metadata,
                         **table_metadata},
            'dataheader': tableheaderstr.split('\t'),
            'data': datastr
        }
    }
    )
    return table_dict

def get_multiplier(datatype, key):
    if key in meas_struct[datatype]['multiplier']:
        return meas_struct[datatype]['multiplier'][key]
    else:
        return 1

def load_tfdata(table_dict):
    obj_list = []
    for f in table_dict:
        datatype = table_dict[f]['meastype']
        if not meas_struct[datatype]:
            raise NotImplementedError
        else:
            for table in table_dict[f]['datatables']:
                t = table_dict[f]['datatables'][table]['data']
                header = table_dict[f]['datatables'][table]['dataheader']
                table_array = np.genfromtxt(t)

                dataobj = meas_struct[datatype]['datatype']()
                dataobj.sample_name = f
                for key, val in meas_struct[datatype]['metadata'].items():
                    v = table_dict[f]['datatables'][table]['metadata'][val]
                    try:
                        m = get_multiplier(datatype, key)
                        dataobj.__setattr__(key, float(v)*m)
                    except TypeError:
                        dataobj.__setattr__(key, v)
                for key, val in meas_struct[datatype]['datatable'].items():
                    m = get_multiplier(datatype, key)
                    # TODO: make this more elegant - E field, current calculated from raw read-in parameters?
                    if key == 'lcm_current': # converts Current density to current
                        dataobj.__setattr__(key, m*table_array[:, val]*dataobj.area)
                    else:
                        dataobj.__setattr__(key, m * table_array[:, val])

                obj_list.append(dataobj)
    return obj_list
