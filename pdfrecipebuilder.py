#!/usr/bin/env python

import collections

import numpy

from pyobjcryst.crystal import Crystal
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe
from diffpy.structure import Lattice

class PDFRecipeBuilder:

    _config_presets = collections.OrderedDict((
        ('rmin', 1.0),
        ('rmax', 20.0),
        ('isotropy', False),
        ('nyquist', False),
        ('fbdelta2', 2.0),
        ('fbbiso', 0.5),
        ('fbqdamp', 0.04),
        ('fbqbroad', 0.03),
    ))


    def __init__(self, crystal, pdfdata, **kw):
        self._config = self._config_presets.copy()
        self.crystal = crystal
        self.pdfdata = dict(pdfdata)
        self.update(**kw)
        return


    def update(self, **kw):
        unknowns = tuple(n for n in kw if not n in self._config)
        if unknowns:
            emsg = ("unsupported keyword argument(s)" +
                    ', '.join(unknowns) + ".")
            raise TypeError(emsg)
        self._config.update(kw)
        return


    def make(self):
        cpdf = PDFContribution('cpdf')
        xobs = self.pdfdata['r']
        yobs = self.pdfdata['g']
        dyobs = self.pdfdata.get('dg')
        cpdf.setObservedProfile(xobs, yobs, dyobs)
        m = self.pdfdata.get('meta', {})
        cpdf.profile.meta.update(m)
        cpdf.addStructure('cif', self.crystal)
        cfg = self._config
        cpdf.setCalculationRange(cfg['rmin'], cfg['rmax'])
        if cfg['nyquist']:
            assert 'qmax' in m, "Nyquist spacing requires 'qmax' value."
            assert m['qmax'] == cpdf.cif.getQmax()
            cpdf.setCalculationRange(dx=numpy.pi / m['qmax'])
        # create FitRecipe
        recipe = FitRecipe()
        recipe.addContribution(cpdf)
        recipe.addVar(cpdf.scale)
        # get symmetry allowed structure parameters
        sgpars = cpdf.cif.phase.sgpars
        # constrain available lattice parameters
        for p in sgpars.latpars:
            recipe.addVar(p, tags=['phase', 'lattice'])
        # constrain free atom positions
        for p in sgpars.xyzpars:
            recipe.addVar(p, tags=['phase', 'positions'])
        # constrain adps
        isosymbol = sgpars.isosymbol
        # get dummy atom with isotropic Biso = fbbiso
        fbbiso = cfg['fbbiso']
        afbiso = _dummyAtomWithUnitBiso(self.crystal)
        afbiso.Bisoequiv = fbbiso
        tags = ['phase', 'adps']
        for p in sgpars.adppars:
            if p.name.startswith(isosymbol):
                recipe.addVar(p, value=p.value or fbbiso, tags=tags)
                continue
            # we have anisotropic site here, but constrain as isotropic
            # if so requested
            if cfg['isotropy']:
                idx = int(p.name.split('_')[-1])
                psc = cpdf.cif.phase.scatterers[idx]
                pbiso = psc.Biso
                n = isosymbol + '_{}'.format(idx)
                v = pbiso.value or fbbiso
                recipe.addVar(pbiso, name=n, value=v, tags=tags)
                continue
            # now we are constraining an anisotropic site here
            # make sure its ADPs are nonzero
            spa = p.par.obj
            if not all((spa.B11, spa.B22, spa.B33)):
                spa.B11 = afbiso.B11
                spa.B22 = afbiso.B22
                spa.B33 = afbiso.B33
                spa.B12 = afbiso.B12
                spa.B13 = afbiso.B13
                spa.B23 = afbiso.B23
            recipe.addVar(p, tags=tags)
        # constrain delta2, qdamp and qbroad
        p = cpdf.cif.delta2
        v = p.value or cfg['fbdelta2']
        recipe.addVar(p, value=v, tag='phase')
        p = cpdf.qdamp
        v = p.value or cfg['fbqdamp']
        recipe.addVar(p, value=v, tag='experiment')
        p = cpdf.qbroad
        v = p.value or cfg['fbqbroad']
        recipe.addVar(p, value=v, tag='experiment')
        return recipe
