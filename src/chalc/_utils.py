import re
from typing import Callable, Any

def interpolate_docstring(names : dict[str, Any] = dict()) -> Callable[[Callable], Callable] :
    def inner(func : Callable) -> Callable:
        docstring : str | None = func.__doc__
        if docstring is not None:
            docstring = docstring.replace('{', '{{').replace('}', '}}')
            docstring = re.sub(r'\${{([^{}]*)}}', r'{\1}', docstring)
            g = { **func.__globals__, **names }
            func.__doc__ = eval('rf"""' + docstring + '"""', g)
        return func
    return inner

