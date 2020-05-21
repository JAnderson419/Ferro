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

def test_tsvload_frequency_parse(input_temp_dhm):
    assert input_temp_dhm.freq == 100