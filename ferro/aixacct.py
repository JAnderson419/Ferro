
import re
import numpy as np
from os.path import basename
from enum import Enum

class MeasEnum(Enum):
    HYSTERESIS = 1
    FATIGUE = 2
    PULSE = 3
    LEAKAGE = 4


def check_datatype(filepath):
    f = open(filepath)
    firstline = f.readline()
    f.close()
    if re.match('DynamicHysteresisResult', firstline):
        return MeasEnum.HYSTERESIS
    elif re.match('Fatigue', firstline):
        return MeasEnum.FATIGUE
    elif re.match('PulseResult', firstline):
        return MeasEnum.PULSE
    elif re.match('LeakageResult'):
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
    table_dict = {}
    tableheaderstr = ''
    table_key = ''

    with open(filepath) as f:
        for line in f:
            line.rstrip('\n')
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
                    if re.match(r'^(/n)?$', line):
                        is_data_table_header = True
                        is_data_table = False
                        table_dict[basename(filepath).split('.')[0]] = {
                            'meastype': datatype,
                            'datatables': {
                                table_key: {
                                    'metadata': {**global_metadata,
                                                 **table_metadata},
                                    'dataheader': tableheaderstr.split('\t'),
                                    'data': datastr
                                }
                            }
                        }
                        table_metadata.clear()
                        datastr = []
                    else:
                        datastr.append(line)
    return table_dict


def load_tfdata(table_dict):

    for f in table_dict:
        print(f)
        for table in table_dict[f]['datatables']:
            t = table_dict[f]['datatables'][table]['data']
            header = table_dict[f]['datatables'][table]['dataheader']
            table_array = np.genfromtxt(t, skip_header=1)
            print(header, table_array)
            exit(0)


if __name__ == '__main__':
    data = read_tfdata(r'C:\Users\Jackson\PycharmProjects\Ferro\tests\testData\RTWhiteB\RTWhiteB_freqs.dat')
    load_tfdata(data)