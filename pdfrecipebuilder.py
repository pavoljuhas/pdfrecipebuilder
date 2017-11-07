#!/usr/bin/env python

from pyobjcryst.crystal import Crystal
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe

class PDFRecipeBuilder:

    def __init__(self, crystal, pdfdata):
        self.crystal = crystal
        self.pdfdata = pdfdata.copy()
        return


    def make(self):
        pass
