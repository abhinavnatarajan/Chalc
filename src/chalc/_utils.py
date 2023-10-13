import re
from typing import Callable

def interpolate_docstring(func : Callable) -> Callable:
    docstring : str | None = func.__doc__
    if docstring is not None:
        docstring = docstring.replace('{', '{{').replace('}', '}}')
        docstring = re.sub(r'\${{([^{}]*)}}', r'{\1}', docstring)
        func.__doc__ = eval('f"""' + docstring + '"""', func.__globals__)
    return func

