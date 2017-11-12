#!/usr/bin/env python

"""Unit tests for pdfrecipefactory.py
"""

import unittest

from pdfrecipefactory import PDFRecipeFactory

# ----------------------------------------------------------------------------

class TestPDFRecipeFactory(unittest.TestCase):

    def setUp(self):
        self.factory = PDFRecipeFactory()
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


    def test_update(self):
        """check PDFRecipeFactory.update()
        """
        return

    def test_make(self):
        """check PDFRecipeFactory.make()
        """
        return

# End of class TestPDFRecipeFactory

# Local Routines -------------------------------------------------------------

def _datafile(filename):
    thisdir = os.path.dirname(os.path.abspath(__file__))
    rv = os.path.join(thisdir, filename)
    return rv


if __name__ == '__main__':
    unittest.main()
