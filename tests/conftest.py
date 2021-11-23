import pytest
import os
import pyanalisi

import faulthandler
faulthandler.enable()


@pytest.fixture(scope='session')
def filepath_tests():
    """Return the absolute filepath of the `tests` folder.

    .. warning:: if this file moves with respect to the `tests` folder, the implementation should change.

    :return: absolute filepath of `tests` folder which is the basepath for all test resources.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(scope='session')
def numpy_traj(filepath_tests):
    import numpy as np
    pos = np.load(filepath_tests + '/data/positions.npy')
    vel = np.load(filepath_tests + '/data/velocities.npy')
    box = np.load(filepath_tests + '/data/cells.npy')
    types = np.zeros(pos.shape[1], dtype = 'i')
    types[-16:-8]=1
    types[-8:]=2
    print('position array shape is {}'.format(pos.shape))
    print('first cell matrix is {}'.format(box[0]))
    return [pos, vel, types, box]

@pytest.fixture(scope='session')
def analisi_traj(numpy_traj):
    return pyanalisi.Trajectory(*numpy_traj,True, False)

@pytest.fixture(scope='session')
def numpy_log(filepath_tests):
    import numpy as np
    with open(filepath_tests + '/data/gk_integral.dat', 'r') as flog:
        headers = flog.readline().split()
    log = np.loadtxt(filepath_tests + '/data/gk_integral.dat', skiprows=1)
    return log, headers

@pytest.fixture(scope='session')
def analisi_log(numpy_log):
    return pyanalisi.ReadLog(*numpy_log)

