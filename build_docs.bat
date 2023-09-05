python docs\stubgen.py chalc --output-dir docs\stubs
xcopy src\chalc\sixpack\ docs\stubs\chalc\sixpack\ /E/Y

sphinx-build -E -b html docs docs\_build