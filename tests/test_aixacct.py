from context import data as hd
from context import aixacct as aix
from os.path import join, dirname, realpath

sampledir = join(dirname(realpath(__file__)), 'testData', 'RTWhiteB')

def test_direct_load_aixACCT_hysteresis_data():
    RTfreqDir = join(sampledir, 'RTWhiteB_freqs')
    RTfreqFiles = hd.dir_read(RTfreqDir)
    RTfreqData = hd.list_read(RTfreqFiles, thickness=255E-7, area=1e-4)

    directreaddict = aix.read_tfdata(join(sampledir, 'RTWhiteB_freqs.dat'))
    directFreqData = aix.load_tfdata(directreaddict)

    assert all(elem in directFreqData for elem in RTfreqData)

def test_direct_load_aixACCT_leakage_data():
    RTlkgDir = join(sampledir, 'RTWhiteB_lkg')
    RTlkgFiles = hd.dir_read(RTlkgDir)
    lkglist = []
    for i in RTlkgFiles:
        l = hd.LeakageData(thickness=255E-7, area=1e-4)
        l.lcm_read(i)
        lkglist.append(l)

    directreaddict = aix.read_tfdata(join(sampledir, 'RTWHITEB_lkg.dat'))
    directlkgData = aix.load_tfdata(directreaddict)
    print(all(elem in directlkgData for elem in lkglist))
    assert all(elem in directlkgData for elem in lkglist)