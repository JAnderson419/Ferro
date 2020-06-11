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

@pytest.fixture
def file_input():
    testdatadir = join(dirname(dirname(realpath(__file__))), "tests", "testData")
    filedir = join(testdatadir, "Typical_AB_Data_RT")
    file = join(filedir, "TypABdata.Hysteresis.2.txt")
    data = hd.HysteresisData()
    data.read_RTHyst(file)
    return data

def test_meta(file_input):
    assert file_input.freq == 100
    assert file_input.thickness == 0.26
    assert file_input.cap_number == 0
    assert file_input.area == 1.00*pow(10,-4)
    assert file_input.temp == 300

def test_data(file_input):
    assert file_input.efield == 346.15
    assert (file_input.voltage[0] == 0.0006) & (file_input.voltage[500] == -0.0006)
    assert (file_input.polarization[0] == -28.230775) & (file_input.polarization[500] == -29.763441)
    assert file_input.vdd == 9.00


@pytest.fixture
def file_lkginput():
    testdatadir = join(dirname(dirname(realpath(__file__))), "tests", "testData")
    filedir = join(testdatadir, "Typical_AB_Data_RT")
    file = join(filedir, "TypABdataLkg.txt")
    data = hd.LeakageData()
    data.read_RTlkg(file)
    return data

def test_meta(file_lkginput):
    assert file_lkginput.thickness == 0.26
    assert file_lkginput.cap_number == 0
    assert file_lkginput.area == 1.00*pow(10,-4)
    assert file_lkginput.temp == 300

def test_data(file_lkginput):
    assert (file_lkginput.lcm_voltage[0] == 9.000244) & (file_lkginput.lcm_voltage[509] == 9.000549)
    assert (file_lkginput.lcm_current[0] == 7.271*pow(10,-10)) & (file_lkginput.lcm_current[509] == 1.596*pow(10,-11))
    assert file_lkginput.vdd == 9.00