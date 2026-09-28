"""Direct Python configuration for the existing Al full-refinement example.

Construction performs no I/O or optimization. Script execution creates and loads
a separate headerless experimental body. Values are independently
transcribed from the reference input and current runner settings, including its
deliberate seven-parameter category selection and non-legacy policies.
"""
from qcbed.configuration import (
    QcbedConfig, MicroscopeConfig, SampleConfig, SpeciesConfig, RegionConfig,
    StoredSliceConfig, AtomConfig, ThermalFactorConfig, ReflectionConfig,
    AbsorptionConfig, AbsorptionSpeciesConfig, ScanConfig, ExperimentalDataConfig,
    RefinementConfig, OptimizerConfig, ParameterBinding, ParameterTarget,
    ParameterCategory as Category, BoundsConfig, TargetBound, CategoryBound,
    RelativeBound, RunConfig, LegacyExportConfig,
)
from qcbed.models import MeshConfig, SliceGeometry
from qcbed.execution import ForwardExecutionConfig


def reference_configuration(*, legacy_header=False) -> QcbedConfig:
    """Default to the headerless copy; legacy_header=True retains old metadata.

    Run this file once to create the separate body file under tmp/. Neither
    constructing nor compiling this configuration performs that conversion.
    """
    reflections = ((1, 1), (1, -1), (-1, 1), (-1, -1),
                   (2, 0), (0, -2), (-2, 0), (2, 0))
    # The final duplicate (2,0) is intentional source data.
    def binding(category, ids, step):
        return ParameterBinding(ParameterTarget(category, ids), (1.05, .95), step)
    parameters = (
        binding(Category.DEPTH, ('specimen',), .5),
        binding(Category.ABSORPTION, ('al', 'b'), .005),
        binding(Category.ABSORPTION, ('al', 'c'), .005),
        binding(Category.DEBYE_WALLER, ('specimen', 'al'), .0002),
        *(binding(Category.SCAN_COORDINATE, (corner, component), .0005)
          for corner in ('A', 'B', 'C') for component in ('first', 'second')),
        *(binding(Category.REFLECTION_MULTIPLIER, ('specimen', f'r{i}', 'elastic'), .0005)
          for i in range(len(reflections))),
    )
    return QcbedConfig(
        title='PACBED EMPAD MONASH E.T. Trondheim',
        microscope=MicroscopeConfig(200., (0, 0, 1)),
        sample=SampleConfig(MeshConfig(16, 16), (SpeciesConfig('al', 'AL', 13, 'AL'),),
            (RegionConfig('specimen', 840.,
                (StoredSliceConfig('al_slice', SliceGeometry(2.0245, 0., 0., 2.8630754, 2.8630754, 90.),
                    (AtomConfig('a0', 'al', 0., 0., .5), AtomConfig('a1', 'al', .5, .5, .5))),),
                (ThermalFactorConfig('al', .86),),
                tuple(ReflectionConfig(f'r{i}', h, k, 1., 1.) for i, (h, k) in enumerate(reflections))),),
            2.5, 0., AbsorptionConfig(True, (AbsorptionSpeciesConfig('al', 0., 1., 1.),), 'YES_ABSORPTION      ')),
        scan=ScanConfig(39, 39, (1., 1., 0.), (1., -1., 0.),
            (-.508156, .407641), (.407641, .508156), (-.407641, -.508156)),
        experiment=ExperimentalDataConfig('in.dat' if legacy_header else 'tmp/configuration_reference/in_body.dat',
            ((0, 0), (1, 1), (1, -1), (-1, -1),
            (-1, 1), (2, 0), (0, -2), (-2, 0), (0, 2)), 1,
            header='validate' if legacy_header else 'absent'),
        refinement=RefinementConfig(
            OptimizerConfig('ANSI', 1, 20, 3, 1e-6, 1e-6, 100., 1., .9, 2026, 1, 1),
            parameters, -1., 'masked', 'reject', (Category.DEPTH, Category.SCAN_COORDINATE),
            BoundsConfig('reject', (TargetBound(ParameterTarget(Category.DEPTH, ('specimen',)), 600., 1200.),),
                tuple(CategoryBound(c, .5, 1.5) for c in
                      (Category.ABSORPTION, Category.DEBYE_WALLER, Category.REFLECTION_MULTIPLIER)),
                (RelativeBound(Category.SCAN_COORDINATE, .2, .05),))),
        execution=ForwardExecutionConfig('batched', 64),
        run=RunConfig('tmp/full_refinement', 'out.dat', 25, 25, 150, 1800.),
        legacy_export=LegacyExportConfig(.5,
            title_raw='PACBED EMPAD MONASH E.T. Trondheim'.ljust(80),
            output_filename_raw='out.dat   ', experimental_filename_raw='in.dat   ' if legacy_header else None),
    )


if __name__ == '__main__':
    from pathlib import Path
    from qcbed import compile_config, load_experimental_data
    from qcbed.legacy_io import _Records

    root = Path(__file__).resolve().parents[1]
    config = compile_config(reference_configuration())
    records = _Records(root/'reference_data/in.dat')
    header = records.read('experimental header')
    actual = tuple(header.integer(name) for name in ('NUMP', 'NUML', 'GRADIENTOPT'))
    expected = (config.state.scan.scan_columns, config.state.scan.scan_rows, config.experiment.gradient_option)
    if actual != expected:
        raise ValueError('reference header differs from canonical configuration')
    destination = root/config.experiment.path
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation protects previous user output; reference data is read-only.
    with destination.open('x', encoding='utf-8') as stream:
        stream.write('\n'.join(records.lines[records.position:])+'\n')
    data = load_experimental_data(config, base_directory=root)
    print(f'{destination}: {data.intensities.shape}, GRADIENTOPT={data.gradient_option}')
