import pytest
from ferro import data as hd
from ferro import aixacct as aix
from os.path import join, dirname, realpath


@pytest.fixture
def input_temp_dhm():
    sampledir = join(dirname(realpath(__file__)), 'testData',  r"hfo2_MFM", "H9_x9y4_1e4_S3_temps")
    file = r'H9 die (9,4) S3 31C 100Hz 3V 1Average Table2.tsv'
    data = hd.HysteresisData()
    data.tsv_read(join(sampledir, file))
    return data


def test_tsvload_temp_parse(input_temp_dhm):
    assert input_temp_dhm.temp == 304



@pytest.fixture
def input_notemp_dhm():
    sampledir = join(dirname(realpath(__file__)), 'testData',  r"hfo2_MFM", "H9_x9y4_1e4_freq")
    file = r'H9 die (9,4) 400Hz 4V 1Average Table18.tsv'
    data = hd.HysteresisData()
    data.tsv_read(join(sampledir, file))
    return data


def test_tsvload_default_temp(input_notemp_dhm):
    assert input_notemp_dhm.temp == 300


def test_tsvload_frequency_parse(input_notemp_dhm):
    assert input_notemp_dhm.freq == 400

def test_hysteresisData_Str(input_notemp_dhm):
    string = input_notemp_dhm.__str__()
    assert string == 'Hysteresis Data, 401 points, -3.97 to 3.96 V, 400.0 Hz, 300K, Pmax = 18.43 uC/cm^2'