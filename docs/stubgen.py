import pybind11_stubgen
from pathlib import Path

if __name__ == '__main__':

    pybind11_stubgen.main()

    stubs_path = (Path(__file__).parents[0].resolve() / 'stubs').resolve()
    stubs = list(stubs_path.rglob('*.pyi'))
    for x in stubs:
        x.replace(x.with_suffix('.py'))
