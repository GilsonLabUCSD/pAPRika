import parmed as pmd

hg = pmd.load_file('tmp.pdb', structure=True)
dum = pmd.topologyobjects.Atom()
dum.name = 'DUM'
dum.mass = 208
dum.atomic_number = 82
dum.xx = -5.0
dum.xy = 0
dum.xz = 0
dum.num = hg.atoms[-1].idx + 1

hg.add_atom(dum, 'DUM', 2)
hg.save('tmp-dum-no-renumber.pdb', renumber=False)


