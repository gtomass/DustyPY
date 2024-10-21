from astroquery.sdss import SDSS
from astropy import coordinates as coords
from astropy.io import fits

# pos = coords.SkyCoord('10h23m19.4847405144s -59d32m04.566918624s', frame='galactic')
# xid = SDSS.query_region(pos, radius='5 arcsec', spectro=True)
pos = coords.SkyCoord('0h8m05.63s +14d50m23.3s', frame='icrs')
xid = SDSS.query_region(pos, radius='5 arcsec', spectro=True)
#print(xid)

sp = SDSS.get_spectra(matches=xid)



print(sp[0][1].header)