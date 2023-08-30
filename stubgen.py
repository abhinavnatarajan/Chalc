import os, pybind11_stubgen
from pybind11_stubgen.structs import(
    Property, 
    QualifiedName, 
    Identifier,
    Docstring,
    Alias,
    Class,
    Method,
    Field)
from pybind11_stubgen.printer import Printer
from pybind11_stubgen.parser.interface import IParser
from pybind11_stubgen.parser.mixins.parse import(
    ExtractSignaturesFromPybind11Docstrings)
from typing import Any
from pathlib import Path
import dataclasses

from pybind11_stubgen.parser.mixins.fix import ReplaceReadWritePropertyWithField

def handle_class_member(
    self, path: QualifiedName, class_: type, obj: Any
) -> Docstring | Alias | Class | list[Method] | Field | Property | None:
    result = super(ReplaceReadWritePropertyWithField, self).handle_class_member(path, class_, obj)    
    return result

ReplaceReadWritePropertyWithField.handle_class_member = handle_class_member

def handle_property(self, path: QualifiedName, prop: Any) -> Property | None:
    result = Property(name=path[-1], modifier=None)

    def get_fake_path(func):
        # Note: pybind *usually* does not include function name
        #       in getter/setter signatures, e.g.:
        #           (arg0: demo._bindings.enum.ConsoleForegroundColor) -> int
        #
        #       Let's pretend the function name is empty if its
        #       docstring starts with `(`
        doc = getattr(func, "__doc__", None)
        if doc is not None and doc.startswith("("):
            return QualifiedName((*path, Identifier("")))
        return path

    if getattr(prop, "__doc__", None) is not None:
        result.doc = prop.__doc__
    if getattr(prop, "fget", None) is not None:
        getters = self.handle_function(get_fake_path(prop.fget), prop.fget)
        if len(getters) == 1:
            result.getter = getters[0]
        if len(getters) > 1:
            raise RuntimeError("Getter overloads")
    if getattr(prop, "fset", None) is not None:
        setters = self.handle_function(get_fake_path(prop.fset), prop.fset)
        if len(setters) == 1:
            result.setter = setters[0]
        if len(setters) > 1:
            raise RuntimeError("Setter overloads")
    if result.getter is None and result.setter is None:
        return None
    return result

ExtractSignaturesFromPybind11Docstrings.handle_property = handle_property

def print_property(self, prop: Property) -> list[str]:
    if not prop.getter:
        # FIXME: support setter-only props
        return []

    # FIXME: add modifier
    result = []

    if prop.doc:
        prop.getter.doc = prop.doc

    result.extend(
        [
            "@property",
            *self.print_function(dataclasses.replace(prop.getter, name=prop.name))
        ]
    )
    
    if prop.setter:
        result.extend(
            [
                f"@{prop.name}.setter",
                *self.print_function(
                    dataclasses.replace(prop.setter, name=prop.name)
                ),
            ]
        )

    return result

Printer.print_property = print_property

if __name__ == '__main__':

    pybind11_stubgen.main()

    # remove existing files to avoid error on windows while renaming
    stubs_path = Path('stubs').resolve()
    old_files = list(stubs_path.rglob('*.py'))
    [os.remove(x) for x in old_files]
    stubs = list(stubs_path.rglob('*.pyi'))
    for x in stubs:
        x.rename(x.with_suffix('.py'))