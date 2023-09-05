import os, pybind11_stubgen
from pathlib import Path

if __name__ == '__main__':

    pybind11_stubgen.main()

    # remove existing files to avoid error on windows while renaming
    stubs_path = (Path(__file__).parents[0].resolve() / 'stubs').resolve()
    stubs = list(stubs_path.rglob('*.pyi'))
    for x in stubs:
        if x.with_suffix('.py').is_file():
            os.remove(x.with_suffix('.py'))
        x.rename(x.with_suffix('.py'))