#!/usr/bin/env python

import numpy

from pyobjcryst.crystal import Crystal
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe
from diffpy.structure import Lattice

class PDFRecipeFactory:

    _config_presets = {
        'rmin' : 1.0,
        'rmax' : 20.0,
        'isotropy' : False,
        'nyquist' : False,
        'fbdelta2' : 2.0,
        'fbbiso' : 0.5,
        'fbqdamp' : 0.04,
        'fbqbroad' : 0.03,
    }


    def __init__(self, **kw):
        for n, v in self._config_presets.items():
            setattr(self, n, v)
        self.update(**kw)
        return


    def update(self, **kw):
        unknowns = tuple(n for n in kw if not n in self._config_presets)
        if unknowns:
            emsg = ("unsupported keyword argument(s)" +
                    ', '.join(unknowns) + ".")
            raise TypeError(emsg)
        for n, v in kw.items():
            setattr(self, n, v)
        return


    def make(self, crystal, r, g, dg=None, meta=None):
        cpdf = PDFContribution('cpdf')
        cpdf.profile.setObservedProfile(r, g, dg)
        m = {} if meta is None else dict(meta)
        cpdf.profile.meta.update(m)
        cpdf.addStructure('cif', crystal)
        cpdf.setCalculationRange(self.rmin, self.rmax)
        if self.nyquist:
            assert 'qmax' in m, "Nyquist spacing requires 'qmax' metadata."
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
        fbbiso = self.fbbiso
        # make a dummy diffpy.structure.Atom with isotropic Biso = fbbiso
        afbiso = _dummyAtomWithBiso(crystal, self.fbbiso)
        tags = ['phase', 'adps']
        for p in sgpars.adppars:
            if p.name.startswith(isosymbol):
                recipe.addVar(p, value=p.value or fbbiso, tags=tags)
                continue
            # we have anisotropic site here, but constrain as isotropic
            # if so requested
            if self.isotropy:
                # extract site index for this sg parameter.  Use it to get
                # the parameter for its Biso value.
                idx = int(p.name.split('_')[-1])
                psite = cpdf.cif.phase.scatterers[idx]
                pbiso = psite.Biso
                n = isosymbol + '_{}'.format(idx)
                v = pbiso.value or fbbiso
                # avoid duplicate constrain
                if recipe.get(n) is None:
                    recipe.addVar(pbiso, name=n, value=v, tags=tags)
                continue
            # here we constrain an anisotropic site.
            # make sure its ADPs are nonzero.
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
        v = p.value or self.fbdelta2
        recipe.addVar(p, value=v, tag='phase')
        p = cpdf.qdamp
        v = p.value or self.fbqdamp
        recipe.addVar(p, value=v, tag='experiment')
        p = cpdf.qbroad
        v = p.value or self.fbqbroad
        recipe.addVar(p, value=v, tag='experiment')
        return recipe

# end of class PDFRecipeFactory


def _dummyAtomWithBiso(crystal, biso):
    from numpy import degrees
    from diffpy.structure import Lattice, Atom
    lattice = Lattice(crystal.a, crystal.b, crystal.c,
                      degrees(crystal.alpha),
                      degrees(crystal.beta),
                      degrees(crystal.gamma))
    a = Atom('X', anisotropy=False, lattice=lattice)
    a.Bisoequiv = biso
    return a
