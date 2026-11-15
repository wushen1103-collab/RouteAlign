"""Connectivity-preserving, symmetry-aware RMSD in the receptor frame.

Only chemical graphs and candidate coordinates are read. All full-graph maps
over resonance forms are enumerated; hitting a bound is an error, not a result.
"""
import numpy as np
from rdkit import Chem

MAX_MAPS = 100000
MAX_RESONANCE = 10000

def heavy(mol):
    m = Chem.RemoveHs(Chem.Mol(mol))
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    if m.GetNumConformers() != 1:
        raise ValueError('exactly one conformer required')
    if not np.isfinite(m.GetConformer().GetPositions()).all():
        raise ValueError('nonfinite candidate coordinates')
    return m

def graph_maps(reference, probe, use_chirality=False):
    if reference.GetNumAtoms() != probe.GetNumAtoms() or reference.GetNumBonds() != probe.GetNumBonds():
        raise ValueError('candidate graph sizes differ')
    resonance = Chem.ResonanceMolSupplier(probe, 0, MAX_RESONANCE + 1)
    nres = len(resonance)
    if nres > MAX_RESONANCE:
        raise RuntimeError('resonance enumeration limit reached')
    maps = set()
    # Sanitize each resonance form to restore aromatic bond flags before
    # matching an aromatic query. The supplier-level matcher alone can omit
    # all benzene maps because its internal forms are Kekulized.
    for form in [probe, *resonance]:
        form = Chem.Mol(form)
        Chem.SanitizeMol(form)
        Chem.AssignStereochemistry(form, cleanIt=True, force=True)
        matches = form.GetSubstructMatches(reference, uniquify=False,
                    useChirality=use_chirality, maxMatches=MAX_MAPS + 1)
        if len(matches) > MAX_MAPS:
            raise RuntimeError('graph-map enumeration limit reached')
        maps.update(tuple(m) for m in matches)
        if len(maps) > MAX_MAPS:
            raise RuntimeError('graph-map enumeration limit reached')
    maps = sorted(maps)
    n = reference.GetNumAtoms()
    if not maps or any(len(m) != n or len(set(m)) != n for m in maps):
        raise ValueError('no complete chemical correspondence')
    # Enforce element and isotope identity explicitly, including resonance maps.
    for mapping in maps:
        for i,j in enumerate(mapping):
            a,b=reference.GetAtomWithIdx(i),probe.GetAtomWithIdx(j)
            if (a.GetAtomicNum(),a.GetIsotope()) != (b.GetAtomicNum(),b.GetIsotope()):
                raise ValueError('element/isotope mismatch')
    return maps,nres

def distance(reference, probe, use_chirality=False):
    a,b=heavy(reference),heavy(probe)
    maps,nres=graph_maps(a,b,use_chirality=use_chirality)
    x=np.asarray(a.GetConformer().GetPositions())
    y=np.asarray(b.GetConformer().GetPositions())
    values=[float(np.sqrt(np.mean(np.sum((x-y[list(m)])**2,axis=1)))) for m in maps]
    best=min(range(len(maps)),key=lambda i:(values[i],maps[i]))
    return values[best],{'mapping':list(maps[best]),'mapping_count':len(maps),
                        'resonance_form_count':nres,'rigid_fit':False,
                        'use_chirality_for_atom_mapping':use_chirality,
                        'isomeric_graph_equal':Chem.MolToSmiles(a)==Chem.MolToSmiles(b)}

def matrix(molecules):
    n=len(molecules);d=np.zeros((n,n));audit=[]
    for i in range(n):
        for j in range(i+1,n):
            v,meta=distance(molecules[i],molecules[j]);d[i,j]=d[j,i]=v
            audit.append({'i':i,'j':j,'rmsd_A':v,**meta})
    return d,audit

def self_test():
    from rdkit.Chem import AllChem
    results={}
    for smiles in ['CC(=O)[O-]','c1ccccc1','C[C@H](O)C(=O)O']:
        m=Chem.AddHs(Chem.MolFromSmiles(smiles));AllChem.EmbedMolecule(m,randomSeed=17);m=heavy(m)
        p=Chem.RenumberAtoms(m,list(reversed(range(m.GetNumAtoms()))))
        v,meta=distance(m,p);assert v<1e-12
        p.GetConformer().SetPositions(np.asarray(p.GetConformer().GetPositions())+[3.,0,0])
        shift,_=distance(m,p);assert abs(shift-3)<1e-10
        results[smiles]={'renumbered_rmsd':v,'translated_rmsd':shift,'maps':meta['mapping_count']}
    # Equivalent terminal oxygens are treated symmetrically across resonance.
    m=Chem.AddHs(Chem.MolFromSmiles('CC(=O)[O-]'));AllChem.EmbedMolecule(m,randomSeed=17);m=heavy(m)
    p=Chem.Mol(m);xyz=np.asarray(p.GetConformer().GetPositions());xyz[[2,3]]=xyz[[3,2]];p.GetConformer().SetPositions(xyz)
    v,meta=distance(m,p);assert v<1e-12 and meta['mapping_count']>=2
    results['resonance_terminal_swap']={'rmsd':v,'maps':meta['mapping_count']}
    m=Chem.AddHs(Chem.MolFromSmiles('C[C@H](O)C(=O)O'));AllChem.EmbedMolecule(m,randomSeed=17);m=heavy(m)
    p=Chem.Mol(m);p.GetAtomWithIdx(1).InvertChirality()
    try: distance(m,p,use_chirality=True)
    except ValueError:results['opposite_stereo_rejected']=True
    else:raise AssertionError('enantiomer incorrectly matched')
    v,audit=distance(m,p)
    assert audit['isomeric_graph_equal'] is False
    results['stereo_tag_drift_recorded_without_coordinate_change']=True
    return results
