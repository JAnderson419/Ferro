from context import data as hd
from context import aixacct as aix
from os.path import join, dirname, realpath


def test_direct_load_aixACCT_hysteresis_data():
    sampledir = join(dirname(realpath(__file__)), 'testData', 'RTWhiteB')
    RTfreqDir = join(sampledir, 'RTWhiteB_Freqs')
    RTfreqFiles = hd.dir_read(RTfreqDir)
    RTfreqData = hd.list_read(RTfreqFiles, thickness=255E-7, area=1e-4)

    directreaddict = aix.read_tfdata(join(sampledir, 'RTWhiteB_freqs.dat'))
    directFreqData = aix.load_tfdata(directreaddict)

    assert all(elem in directFreqData for elem in RTfreqData)
