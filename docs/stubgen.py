import os, pybind11_stubgen
from pathlib import Path

if __name__ == '__main__':

    pybind11_stubgen.main()

    # remove existing files to avoid error on windows while renaming
    stubs_path = (Path(__file__).parents[0].resolve() / 'stubs').resolve()
    old_files = list(stubs_path.rglob('*.py'))
    [os.remove(x) for x in old_files]
    stubs = list(stubs_path.rglob('*.pyi'))
    for x in stubs:
        x.rename(x.with_suffix('.py'))