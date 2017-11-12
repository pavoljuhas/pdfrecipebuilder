#!/usr/bin/env python

"""Unit tests for pdfrecipefactory.py
"""

import os
import unittest

from pyobjcryst import loadCrystal
from diffpy.structure import loadStructure
from diffpy.utils.parsers import loadData
from pdfrecipefactory import PDFRecipeFactory

# ----------------------------------------------------------------------------

class TestPDFRecipeFactory(unittest.TestCase):

    def setUp(self):
        self.factory = PDFRecipeFactory()
        self.cifbto = _datafile('BaTiO3-Pm3m.cif')
        data = loadData(_datafile('BaTiO3.gr'))
        self.r, self.g = data.T
        return


    def test___init__(self):
        """check PDFRecipeFactory.__init__()
        """
        self.assertEqual(0.04, self.factory.fbqdamp)
        self.assertEqual(20, self.factory.rmax)
        rf2 = PDFRecipeFactory(rmax=13, isotropy=True)
        self.assertEqual(13, rf2.rmax)
        self.assertTrue(rf2.isotropy)
        self.assertRaises(TypeError, PDFRecipeFactory, xinvalid=4)
        return


    def test_make_std(self):
        """check PDFRecipeFactory.make() in standard setup.
        """
        crst = loadCrystal(self.cifbto)
        meta = dict(qmax=26, stype='N')
        recipe = self.factory.make(crst, self.r, self.g, meta=meta)
        self.assertEqual(1, recipe.cpdf.r.value[0])
        self.assertEqual(20, recipe.cpdf.r.value[-1])
        # check the effect of meta entries
        calc = recipe.cpdf.cif._calc
        self.assertEqual(26, calc.qmax)
        self.assertEqual('N', calc.scatteringfactortable.radiationType())
        # should have 4 parameters for B factors, 2 are isotropic
        bnames = [n for n in recipe.names if n.startswith('B')]
        self.assertEqual(4, len(bnames))
        self.assertEqual(2, len([n for n in bnames if n.startswith('Biso')]))
        # test type checking
        stru = loadStructure(self.cifbto)
        self.assertRaises(TypeError, self.factory.make,
                          stru, self.r, self.g, meta)
        return


    def test_make_biso_fallback(self):
        """check if Biso fallback in PDFRecipeFactory is applied.
        """
        self.factory.fbbiso = 0.37
        crst = loadCrystal(self.cifbto)
        # zero all B values in the model
        spreg = crst.GetScatteringPowerRegistry()
        for i in range(spreg.GetNb()):
            sp = spreg.GetObj(i)
            sp.B11 = sp.B22 = sp.B33 = 0.0
        recipe = self.factory.make(crst, self.r, self.g)
        spba = spreg.GetObj(0)
        fbbiso = self.factory.fbbiso
        self.assertEqual(fbbiso, spba.Biso)
        spo = spreg.GetObj(2)
        self.assertEqual(fbbiso, spo.B11)
        self.assertEqual(fbbiso, spo.B22)
        self.assertEqual(fbbiso, spo.B33)
        return


    def test_make_isotropy(self):
        """check PDFRecipeFactory.make() with `isotropy` option.
        """
        crst = loadCrystal(self.cifbto)
        self.factory.isotropy = True
        recipe = self.factory.make(crst, self.r, self.g)
        bnames = [n for n in recipe.names if n.startswith('B')]
        self.assertEqual(3, len(bnames))
        self.assertTrue(all(n.startswith('Biso') for n in bnames))
        self.assertEqual(self.factory.fbbiso, recipe.Biso_0.value)
        return


    def test_make_nyquist(self):
        """check PDFRecipeFactory.make() with `nyquist` option.
        """
        crst = loadCrystal(self.cifbto)
        self.factory.nyquist = True
        self.assertRaises(Exception, self.factory.make,
                          crst, self.r, self.g)
        meta = dict(qmax=26)
        recnq = self.factory.make(crst, self.r, self.g, meta=meta)
        rgrid = recnq.cpdf.r.value
        self.assertTrue(rgrid[1] - rgrid[0] > 0.1)
        return

# End of class TestPDFRecipeFactory

# Local Routines -------------------------------------------------------------

def _datafile(filename):
    thisdir = os.path.dirname(os.path.abspath(__file__))
    rv = os.path.join(thisdir, filename)
    return rv


if __name__ == '__main__':
    unittest.main()
