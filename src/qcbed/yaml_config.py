"""Optional, strict YAML frontend; the typed model and compiler own semantics."""
from dataclasses import MISSING, fields, is_dataclass, replace
from enum import Enum
from pathlib import Path
from types import UnionType
from typing import get_args, get_origin, get_type_hints

from .configuration import QcbedConfig, _typed
from .config_compiler import compile_config


def _decode(value, annotation, path):
    origin, args = get_origin(annotation), get_args(annotation)
    if origin is UnionType:
        if value is not None and type(None) in args:
            return _decode(value,next(t for t in args if t is not type(None)),path)
        for option in args:
            try:
                return _decode(value, option, path)
            except (TypeError, ValueError):
                pass
        raise ValueError(f'{path}: expected {annotation}')
    if origin is tuple:
        if not isinstance(value, list):
            raise ValueError(f'{path}: expected a YAML sequence')
        types = (args[0],)*len(value) if len(args) == 2 and args[1] is Ellipsis else args
        if len(types) == 3 and len(value) == 4 and ('hkl' in path or 'zone_axis' in path):
            raise ValueError(f'{path}: 4-index Miller–Bravais notation is not supported; provide canonical 3-index input')
        if len(types) != len(value):
            raise ValueError(f'{path}: expected {len(types)} entries')
        return tuple(_decode(v,t,f'{path}[{i}]') for i,(v,t) in enumerate(zip(value,types)))
    if is_dataclass(annotation):
        if not isinstance(value, dict):
            raise ValueError(f'{path}: expected a mapping')
        declarations = {f.name:f for f in fields(annotation)}
        unknown = set(value)-set(declarations)
        if unknown:
            raise ValueError(f'{path}: unknown fields {sorted(unknown)}')
        for name, field in declarations.items():
            if name not in value and field.default is MISSING and field.default_factory is MISSING:
                raise ValueError(f'{path}.{name}: missing required field')
        hints = get_type_hints(annotation)
        return annotation(**{k:_decode(v,hints[k],f'{path}.{k}') for k,v in value.items()})
    if isinstance(annotation,type) and issubclass(annotation,Enum):
        for member in annotation:
            if type(value) is type(member.value) and value == member.value:
                return member
        raise ValueError(f'{path}: invalid {annotation.__name__} value {value!r}')
    return _typed(value,annotation,path)


def read_yaml_document(path):
    """Strict inert document parsing shared with an unprepared GUI draft."""
    import yaml

    class Loader(yaml.SafeLoader):
        def compose_node(self, parent, index):
            if self.check_event(yaml.AliasEvent):
                raise ValueError('YAML aliases are not supported')
            return super().compose_node(parent,index)

    def mapping(loader,node):
        result = {}
        for key_node,value_node in node.value:
            key = loader.construct_object(key_node,deep=True)
            if type(key) is not str:
                raise ValueError(f'{key_node.start_mark}: mapping keys must be strings')
            if key in result:
                raise ValueError(f'{key_node.start_mark}: duplicate YAML key {key!r}')
            result[key] = loader.construct_object(value_node,deep=True)
        return result

    Loader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,mapping)
    source = Path(path).resolve()
    with source.open(encoding='utf-8') as stream:
        value=yaml.load(stream,Loader=Loader)
    return value


def load_yaml_config(path: str | Path) -> QcbedConfig:
    """Load one strict YAML document; resolve paths against its directory.

    PyYAML is imported only on use. Safe scalar/mapping/sequence tags only;
    aliases and merge keys are rejected. No expressions or substitutions.
    Legacy compilation remains I/O-free. Opt-in crystallography validates its
    external structure through the explicit source-loading frontend.
    """
    value=read_yaml_document(path)
    source=Path(path).resolve()
    if isinstance(value,dict) and any(k in value for k in ('data_prep','atomic_model','qcbed')):
        from .modern_yaml import decode_modern
        config=decode_modern(value,source.parent)
    else:config=_decode(value,QcbedConfig,'config')
    from .config_paths import resolve_config_paths
    resolved = resolve_config_paths(config, base_directory=source.parent)
    if resolved.crystallography is not None:
        from .crystallography_config import compile_config_from_sources
        compile_config_from_sources(resolved)
    return resolved
