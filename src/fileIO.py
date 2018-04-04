### CLASSIC_EBM
### Jake Aylmer
###
### Functions for file input/output and saving figures/data.
### ---------------------------------------------------------------------------

from __future__ import division
import sys, os, numpy as np, matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

def SaveFigures(figures, subdir, ext='.pdf'):
    """Save several figures to plots\subdir with filenames given by their
    canvas window titles. Makes the directory if it does not exist and warns/
    aborts if about to over-write a file.
    
    --Args--
    figures : list of MatPlotLib figure objects.
    subdir  : string, name of sub-directory to which figures should be saved.
    (ext)   : string, type of file to save as its extension (default '.pdf').
    """
    dirname = os.path.join(os.path.dirname(__file__), '..', 'plots', subdir)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)
    passall = False
    for f in figures:
        filename_to_save = str(f.canvas.manager.window.windowTitle()) + ext
        if passall:
            f.savefig(os.path.join(dirname, filename_to_save))
        else:
            if os.path.isfile(os.path.join(dirname, filename_to_save)):
                prt='About to overwrite \"%s\"! Continue? [Y / Y-ALL / N] > '%(
                    filename_to_save)
                check = raw_input(prt)
                if not (check == 'Y' or check == 'Y-ALL'):
                    pass
                else:
                    f.savefig(os.path.join(dirname, filename_to_save))
                passall = check=='Y-ALL'
            else:
                f.savefig(os.path.join(dirname, filename_to_save))
    pass
