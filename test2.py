
from astropy.table import Table

if __name__ == "__main__":
    radius=5
    target="HD302821"
    #target="HD12"
    sed=Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={target}&-c.rs={radius}")
    #http://vizier.cds.unistra.fr/vizier/sed/?submitSimbad=Photometry&-c=10+23+19.4847405144-59+32+04.566918624&-c.r=5&-c.u=arcsec&show_settings=1
    #sed = Table.read(f'http://vizier.cds.unistra.fr/vizier/sed/?submitSimbad=Photometry&-c={target}&-c.r={radius}&-c.u=arcsec', format='csv')
    #sed= Table.read(f"https://vizier.cds.unistra.fr/viz-bin/sed?-c={target}&-c.rs={radius}")
    print(sed)


    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure()
    ax = fig.add_subplot()
    ax.set_xlabel('freq (GHz)')
    ax.set_ylabel('Flux (Jsky)')
    
    ax.scatter(sed["sed_freq"], sed["sed_flux"])
    plt.show()
