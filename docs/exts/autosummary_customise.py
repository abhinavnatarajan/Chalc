__version__ = '0.1.0'
from sphinx.ext.autosummary.generate import(
    AutosummaryRenderer, 
    ModuleScanner, 
    _get_members, 
    _get_module_attrs,
    _get_modules)
from sphinx.ext.autosummary import get_documenter
from typing import Any

# autosummary does not correctly get parse class properties
# this fixes that problem
def generate_autosummary_content(name: str, obj: Any, parent: Any,
                                 template: AutosummaryRenderer, template_name: str,
                                 imported_members: bool, app: Any,
                                 recursive: bool, context: dict,
                                 modname: str | None = None,
                                 qualname: str | None = None) -> str:
    doc = get_documenter(app, obj, parent)

    ns: dict[str, Any] = {}
    ns.update(context)

    if doc.objtype == 'module':
        scanner = ModuleScanner(app, obj)
        ns['members'] = scanner.scan(imported_members)

        respect_module_all = not app.config.autosummary_ignore_module_all
        imported_members = imported_members or ('__all__' in dir(obj) and respect_module_all)

        ns['functions'], ns['all_functions'] = \
            _get_members(doc, app, obj, {'function'}, imported=imported_members)
        ns['classes'], ns['all_classes'] = \
            _get_members(doc, app, obj, {'class'}, imported=imported_members)
        ns['exceptions'], ns['all_exceptions'] = \
            _get_members(doc, app, obj, {'exception'}, imported=imported_members)
        ns['attributes'], ns['all_attributes'] = \
            _get_module_attrs(name, ns['members'])
        ispackage = hasattr(obj, '__path__')
        if ispackage and recursive:
            # Use members that are not modules as skip list, because it would then mean
            # that module was overwritten in the package namespace
            skip = (
                ns["all_functions"]
                + ns["all_classes"]
                + ns["all_exceptions"]
                + ns["all_attributes"]
            )

            # If respect_module_all and module has a __all__ attribute, first get
            # modules that were explicitly imported. Next, find the rest with the
            # get_modules method, but only put in "public" modules that are in the
            # __all__ list
            #
            # Otherwise, use get_modules method normally
            if respect_module_all and '__all__' in dir(obj):
                imported_modules, all_imported_modules = \
                    _get_members(doc, app, obj, {'module'}, imported=True)
                skip += all_imported_modules
                imported_modules = [name + '.' + modname for modname in imported_modules]
                all_imported_modules = \
                    [name + '.' + modname for modname in all_imported_modules]
                public_members = getall(obj)
            else:
                imported_modules, all_imported_modules = [], []
                public_members = None

            modules, all_modules = _get_modules(obj, skip=skip, name=name,
                                                public_members=public_members)
            ns['modules'] = imported_modules + modules
            ns["all_modules"] = all_imported_modules + all_modules
    elif doc.objtype == 'class':
        ns['members'] = dir(obj)
        ns['inherited_members'] = \
            set(dir(obj)) - set(obj.__dict__.keys())
        ns['methods'], ns['all_methods'] = \
            _get_members(doc, app, obj, {'method'}, include_public={'__init__'})
        ns['attributes'], ns['all_attributes'] = \
            _get_members(doc, app, obj, {'attribute'})
        ns['properties'], ns['all_properties'] = \
            _get_members(doc, app, obj, {'property'})
        ns['classes'], ns['all_classes'] = \
            _get_members(doc, app, obj, {'class'})

    if modname is None or qualname is None:
        modname, qualname = _split_full_qualified_name(name)

    if doc.objtype in ('method', 'attribute', 'property'):
        ns['class'] = qualname.rsplit(".", 1)[0]

    if doc.objtype in ('class',):
        shortname = qualname
    else:
        shortname = qualname.rsplit(".", 1)[-1]

    ns['fullname'] = name
    ns['module'] = modname
    ns['objname'] = qualname
    ns['name'] = shortname

    ns['objtype'] = doc.objtype
    ns['underline'] = len(name) * '='

    if template_name:
        return template.render(template_name, ns)
    else:
        return template.render(doc.objtype, ns)


def autosummary_modify_content(func):
    import sphinx.ext.autosummary

    sphinx.ext.autosummary.generate.generate_autosummary_content = func

def setup(app):
    app.setup_extension("sphinx.ext.autosummary")
    autosummary_modify_content(generate_autosummary_content)

    return {"version": __version__, "parallel_read_safe": True}