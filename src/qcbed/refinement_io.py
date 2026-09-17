"""Refinement records layered over the existing FOR007 prefix reader."""

from dataclasses import replace
from pathlib import Path

from .legacy_io import _Records, _parse_for007_records
from .scan import parse_for007_scan
from .refinement import (RefinementRegion, RefinementReflection, RefinementState,
                         build_refinement_parameter_layout)


class _CaptureRead:
    def __init__(self, read, values):
        self.read, self.values = read, values

    def __getattr__(self, name):
        method = getattr(self.read, name)
        if name not in ("real", "integer", "flag"):
            return method

        def capture(*args, **kwargs):
            value = method(*args, **kwargs)
            self.values.append(value)
            return value
        return capture


class _CaptureRecords(_Records):
    """Retain validated prefix fields without duplicating its READ syntax."""
    def __init__(self, path):
        super().__init__(path)
        self.values = {}

    def read(self, section):
        values = self.values.setdefault(section, [])
        return _CaptureRead(super().read(section), values)


def parse_refinement_state(path, *, experimental_header_path=None) -> RefinementState:
    """Read the refinement state; preserve the existing controls-consuming API."""
    return _parse_refinement_input(path, experimental_header_path=experimental_header_path)[0]


def _parse_refinement_input(path, *, experimental_header_path=None):
    """Read active masks, physical values and D1/D2 without running refinement.

    Existing FOR007 and scan parser behavior is unchanged. Optimizer controls
    and PLAMDA records are consumed only to locate the range records.
    """
    records = _CaptureRecords(Path(path))
    model = _parse_for007_records(records)
    scan = parse_for007_scan(path, experimental_header_path=experimental_header_path)
    v = records.values
    depths = v["region depths and masks"]
    regions = []
    for i, count in enumerate(v["region slice counts"]):
        flag = v[f"region {i+1} vacancy flag"][0]
        vacancy = v[f"region {i+1} vacancy parameters"]
        regions.append(RefinementRegion(depths[2*i], depths[2*i+1], count, flag,
            tuple(vacancy[::2]), tuple(vacancy[1::2]),
            tuple(v[f"Debye-Waller factors region {i+1}"]),
            tuple(v[f"Debye-Waller masks region {i+1}"]),
            tuple(v[f"free-surface parameters region {i+1}"])))
    absorption, masks = [], []
    if model.absorption_enabled:
        for i in range(len(model.species)):
            row = v[f"absorption coefficients {i+1}"]
            absorption.append(tuple(row[::2]))
            masks.append(tuple(row[1::2]))
    read = records.read("refinement scan basis")
    for _ in range(6):
        read.real("scan basis component")
    for name in ("A", "B", "C"):
        read = records.read(f"refinement corner {name}")
        for _ in range(2):
            read.real("scan corner")

    def mask(section):
        raw = records.next(section)
        value = raw[:2].ljust(2)
        if value not in ("NN", "YN", "NY", "YY"):
            raise records.error(section, raw, "NN, YN, NY or YY in FORMAT(A2)")
        return value

    scan_masks = tuple(mask(f"scan mask {i+1}") for i in range(3))
    iref = records.read("IREF").real("IREF")
    reflections, symmetry = [], []
    for i in range(len(regions)):
        count = records.read(f"region {i+1} reflection count").integer("NRF", minimum=0)
        raw = records.next(f"region {i+1} symmetry")
        sym = raw[:4].ljust(4)
        if sym not in ("NONE", "INVE", "SQUA", "RECT", "HEXA"):
            raise records.error("symmetry", raw, "NONE, INVE, SQUA, RECT or HEXA")
        symmetry.append(sym)
        row = []
        for j in range(count):
            read = records.read(f"region {i+1} reflection {j+1}")
            row.append((read.integer("H"), read.integer("K"), read.real("SFFAC"), read.real("ABSFAC")))
        reflections.append(tuple(RefinementReflection(*item, mask(f"region {i+1} reflection mask {j+1}"))
                                 for j, item in enumerate(row)))
    state = RefinementState(model, tuple(regions), *v["lattice dilation"], tuple(absorption), tuple(masks),
        scan, scan_masks, tuple(reflections), tuple(symmetry), v["DSAP"][0], iref)
    n = len(build_refinement_parameter_layout(state).parameters)
    read = records.read("optimizer selector (not executed)")
    controls = {"method": read.token().value, "paraopt": read.integer("PARAOPT"),
                "loud": read.integer("LOUD")}
    control_record = records.position
    for name in ("LOOPS", "MAXFSIMP", "MAXFPARA"):
        controls[name.lower()] = records.read(name).integer(name)
    for name in ("FTOL", "CHIMIN", "MAXTEMP", "MINTEMP", "COOLRATE"):
        controls[name.lower()] = records.read(name).real(name)
    # Empty implied-DO READ still consumes a record in the source.
    plamda = []
    if n:
        read = records.read("PLAMDA (not applied)")
        for _ in range(n):
            plamda.append(read.real("PLAMDA"))
    else:
        records.next("empty PLAMDA record")
    ranges = []
    for i in range(n):
        read = records.read(f"parameter {i+1} range factors")
        ranges.append((read.real("D1"), read.real("D2")))
    controls["plamda"] = tuple(plamda)
    return replace(state, ranges=tuple(ranges)), controls, records, control_record
