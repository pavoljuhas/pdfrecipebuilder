from diffpy.srfit.fitbase import FitRecipe
from diffpy.srfit.pdf import PDFContribution
from pyobjcryst import loadCrystal
from diffpy.structure import loadStructure

config = {
    'rmin' : 1,
    'rmax' : 20,
    'Biso_default' : 0.5,
    'delta2_default' : 2.0,
}

crst = loadCrystal('BaTiO3.cif')

cpdf = PDFContribution('cpdf')
cpdf.loadData('BaTiO3.gr')
cpdf.addStructure('cif', crst)
cpdf.profile.setCalculationRange(1, 20)

recipe = FitRecipe()
recipe.addContribution(cpdf)
recipe.addVar(cpdf.scale)
for p in cpdf.cif.phase.sgpars.latpars:
    recipe.addVar(p, tags=['phase', 'lattice'])
for p in cpdf.cif.phase.sgpars.xyzpars:
    recipe.addVar(p, tags=['phase', 'positions'])
for p in cpdf.cif.phase.sgpars.adppars:
    if p.value == 0:
        spa = p.par.obj
        p.value = spa.Biso or spa.B11 or config['Biso_default']
    recipe.addVar(p, tags=['phase', 'adps'])
recipe.addVar(cpdf.cif.delta2, value=config['delta2_default'], tag='phase')
recipe.addVar(cpdf.qdamp, value=cpdf.qdamp.value or 0.04, tag='experiment')
