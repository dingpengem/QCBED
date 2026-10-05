"""Real CIF thermal ownership; no objective or optimizer evaluations."""
from dataclasses import replace
from pathlib import Path
import pytest
from qcbed import compile_config_from_sources
from qcbed.configuration import ThermalFactorConfig
from test_modern_grouped_run import modern_config


def test_real_cif_site_group_b_needs_no_legacy_duplicate(tmp_path):
    c = modern_config(tmp_path)
    c = replace(c, crystallography=replace(c.crystallography,
        structure=replace(c.crystallography.structure,
            path=str(Path(__file__).resolve().parents[1]/'reference_data/Al_FCC.cif'),format=None)))
    assert c.crystallography.debye_waller == ()
    grouped = compile_config_from_sources(c)
    assert grouped.crystallography.site_groups[0].id == 'Al_fcc'
    b = c.crystallography.site_groups[0].debye_waller_b
    assert all(s.debye_waller_b == b for s in grouped.crystallography.original_structure.sites)
    legacy = replace(c, crystallography=replace(c.crystallography,
        site_groups=(),debye_waller=(ThermalFactorConfig('Al', b),)),
        refinement=replace(c.refinement,parameters=(),
            bounds=replace(c.refinement.bounds,specific=())))
    assert compile_config_from_sources(legacy).state.regions[0].debye_waller == grouped.state.regions[0].debye_waller
    absent = replace(legacy,crystallography=replace(legacy.crystallography,debye_waller=()))
    with pytest.raises(ValueError,match='Site-group configs do not require'):
        compile_config_from_sources(absent)
