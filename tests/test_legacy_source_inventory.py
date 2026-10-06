"""Current readable-source integrity, not compiled executable equivalence."""
import importlib.util
import json
from pathlib import Path
from private_reference import private_reference_root

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("source_inventory", ROOT/"tools/legacy_source_inventory.py")
inventory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inventory)


def test_clean_sources_and_unique_project_definitions():
    m = inventory.build_manifest(private_reference_root())
    assert len(m['files']) == 14
    assert all(f['readable'] and not f['nul'] and not f['truncation_flags'] for f in m['files'])
    assert all(len(v) == 1 for v in m['routines'].values())
    assert not [c for c in m['calls'] if c['resolution'] == 'definition missing']
    assert not [c for c in m['calls'] if c['signature_status'] == 'argument-count mismatch']
    assert [Path(f['file']).name for f in m['files'] if not f['units']] == ['MULTIS.FTN']


def test_critical_signatures_and_checked_manifest():
    m = inventory.build_manifest(private_reference_root())
    expected = {
        'ATOMICSCATFACT':'ATNUM,STL,OPTION', 'FFT':'ISIGN',
        'CDFT':'N,ISGN,A,IP,W', 'CDFT2D':'N1MAX,N1,N2,ISGN,A,T,IP,W',
        'ONEIT':'PT','ONEITPARA':'PT','GRADIENT':'','CHISQCAL':'CHI2',
        'ARRANGEPARS':'','UPDPARS':'PT','ANNEAL':'GOF,MAXITS',
        'AMOEBA':'P,Y,MP,NP,NDIM,CHISQ,FTOL,MXF,FMIN,NUMF,SCAL',
        'PARAMIN':'BOT,TOP,GOF,TOLER',
        'SIMPARA':'GOF,MAXITS,HITEMP,LOTEMP,COOLING,UNCERT,INVAR,OUTVAR,NPARAM,POPT,NFUNC,VERB,OUNIT',
        'RANDO':'POW','NRANDO':'GRAIN,CENTRAL,SIG','URANDO':'GRAIN,CENTRAL,SIG','UPDATE':''}
    for name,args in expected.items():
        assert m['routines'][name][0]['arguments'] == (args.split(',') if args else [])
    assert 'FFT2D' not in m['routines']  # actual 2-D entry is CDFT2D
    assert 'SCAL' not in m['routines']   # shared parameter-coordinate array
    saved = json.loads((private_reference_root()/'docs/legacy_source_manifest.json').read_text(encoding='utf-8'))
    semantic = lambda x: {k:[(u['file'],u['type'],u['arguments']) for u in v] for k,v in x['routines'].items()}
    assert semantic(m) == semantic(saved)  # harmless line/whitespace changes allowed


def test_continuations_nested_arguments_and_comment_exclusion():
    text = "C SUBROUTINE FAKE(A)\n      subroutine demo(a,\n     * b)\n\tCALL worker(a(1,2), DBLE(b)) ! comment\n      END\n"
    s = list(inventory.statements(text))
    assert inventory.UNIT.match(s[0]['text'])[3].upper() == 'DEMO'
    assert s[0]['line'] == 2 and s[0]['end_line'] == 3
    assert inventory.split_args("a(1,2), DBLE(b), 'a,b'") == ['a(1,2)', 'DBLE(b)', "'a,b'"]
    assert s[1]['text'] == 'CALL worker(a(1,2), DBLE(b))'


def test_known_abi_hazards_are_visible_not_silently_declared_safe():
    m = inventory.build_manifest(private_reference_root())
    assert any(c['caller']=='ANNEAL' and c['callee']=='NRANDO' and c['visible_kind_disagreements'] for c in m['calls'])
    assert any(c['caller']=='MULTIS' and c['callee']=='AMOEBA' and c['visible_kind_disagreements'] for c in m['calls'])
    assert any(c['caller']=='SIMPARA' and c['callee']=='NRANDO' and c['visible_kind_disagreements'] for c in m['calls'])
