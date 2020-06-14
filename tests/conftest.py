import pytest
import os

@pytest.fixture(scope='session')
def filepath_tests():
    """Return the absolute filepath of the `tests` folder.

    .. warning:: if this file moves with respect to the `tests` folder, the implementation should change.

    :return: absolute filepath of `tests` folder which is the basepath for all test resources.
    """
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture
def numpy_traj(filepath_tests):
    import numpy as np
    pos = np.load(filepath_tests + '/data/positions.npy')
    vel = np.load(filepath_tests + '/data/velocities.npy')
    box = np.load(filepath_tests + '/data/cells.npy')
    types = np.zeros(pos.shape[1], dtype = np.int32)
    types[-8:]=1
    return [pos, vel, types, box]

@pytest.fixture
def analisi_traj(numpy_traj):
    import pyanalisi
    return pyanalisi.Trajectory(*numpy_traj,True, False)

