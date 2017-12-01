#!/usr/bin/env python

import numpy

from pyobjcryst.crystal import Crystal
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe

# ----------------------------------------------------------------------------

class PDFRecipeFactory:
    """
    Factory that makes srfit recipe from CIF structure and observed PDF.

    The CIF structure is passed as a pyobjcryst.Crystal object.

    Attributes
    ----------
    rmin : float
        The lower r-bound for the fitted PDF range.
    rmax : float
        The upper r-bound for the fitted PDF range.
    isotropy : bool
        Flag for using isotropic displacement parameters at all sites
        even when the space group allows for anisotropic displacement.
        The default is ``False``.
    nyquist : bool
        Flag for fitting PDF data re-sampled to Nyquist r-spacing.
        When set the ``qmax`` value must be passed in the `make`
        method `meta` argument.  The default is False.
    fbdelta2 : float
        Fallback value for the *delta2* sharpening factor.
    fbbiso : float
        Fallback value for the site isotropic displacement parameters
        Biso when missing or zero in the CIF structure.
    fbqdamp : float
        Fallback value for the q-resolution damping factor *qdamp*.
    fbqbroad : float
        Fallback value for the resolution-related peak broadening.

    Note
    ----
    The initial `delta2`, `qdamp`, `qbroad` variables can be set via
    the `meta` argument of the `make` method.  Otherwise the recipe
    starts with their fallback values.

    Parameters
    ----------
    kw : keywords
        The keyword arguments to set the factory configuration.
        These can be the instance attributes listed above.
    """

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
        """
        Parameters
        ----------
        kw : keywords
            The keyword arguments to update the factory configuration.
            These can be the attributes listed in the class docstring.
        """
        unknowns = tuple(n for n in kw if not n in self._config_presets)
        if unknowns:
            emsg = ("unsupported keyword argument(s)" +
                    ', '.join(unknowns) + ".")
            raise TypeError(emsg)
        for n, v in kw.items():
            setattr(self, n, v)
        return


    def make(self, crystal, r, g, dg=None, meta=None):
        """
        Construct new srfit recipe from CIF structure and PDF data

        Parameters
        ----------
        crystal : pyobjcryst.Crystal
            The CIF structure to be fitted to the PDF data in-place.
        r : array_like
            The r-grid of the fitted PDF dataset in Angstroms.
        g : array_like
            The fitted PDF values per each `r` point.
        dg : array_like, optional
            The estimated standard deviations at each of `g` values.
            When unspecified, *dg* is assumed 1 leading to underestimated
            standard errors of the refined variables.
        meta : dict, optional
            A dictionary of extra metadata to be used when constructing
            `PDFContribution` in the srfit recipe.  The common recognized
            keys are ``stype`` for radiation type ("X" or "N"), ``qmax``
            for the Q-range used in the experiment, ``delta1``, ``delta2``,
            ``qbroad`` for peak sharpening and broadening factors and
            ``qdamp`` for the Q-resolution related signal dampening.

        Returns
        -------
        recipe : FitRecipe
            The new FitRecipe for in-place fitting of `crystal` to PDF data.
        """
        if not isinstance(crystal, Crystal):
            emsg = "crystal must be of the pyobjcryst.Crystal type."
            raise TypeError(emsg)
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
        recipe.clearFitHooks()
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
                # avoid applying duplicate constrain to pbiso
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

# Local Routines -------------------------------------------------------------

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
