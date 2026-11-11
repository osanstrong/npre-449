import sys
args = sys.argv
if len(args) < 3:
    raise ValueError("Error: missing arguments (python scrubPython.py <input> <output>)")
# fpath = "./ostrong2_NPRE449_cp.py"
# outpath = "./report/cp_export.py"
fpath = args[1]
outpath = args[2]

with open(fpath, 'r') as scriptfile:
    text = scriptfile.read().replace('²', '2').replace('ρ', 'rho').replace('χ', 'chi').replace('⁰', 'deg').replace('μ', 'mu').replace('π', 'pi')
    with open(outpath, 'w') as outfile:
        outfile.write(text)