import pytest
import numpy as np
from testbook import testbook


@pytest.fixture(scope='module')
def tb():
    with testbook('../notebooks/calc_inspector.ipynb', execute=True) as tb:
        yield tb

def test_density(tb,num_regression):
    tb.inject("_cell=density[1].flatten().tolist()")
    tb.inject("_dens=density[0].flatten().tolist()")
    density = tb.ref("_dens")
    cell = tb.ref("_cell")
    num_regression.check({
          'cell':np.array(density,dtype=float),
          'density':np.array(cell)})
    

def disabledtest_sh(tb,num_regression):
    tb.inject("_shfit=np.array(shp[3]).tolist()")
    shfit = tb.ref("_shfit")
    num_regression.check({
          'shfit':np.array(shfit).flatten()
             })

def test_gofr(tb,num_regression):
    tb.inject("_gofr=gofr.tolist()")
    gofr=tb.ref("_gofr")
    num_regression.check({'gofr':np.array(gofr).flatten()})
