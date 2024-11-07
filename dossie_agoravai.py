from astropy.table import Table
import os
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from pycasso2 import FitsCube
from pycasso2.geometry import find_peak

plt.rcParams['image.origin'] = 'lower'

def draw_r50(ax, x0, y0, a, b_a, PA, color='limegreen', lw=1.5):
    from matplotlib.patches import Ellipse
    center = np.array([x0, y0])
    theta = PA - 90
    e = Ellipse(center, height=2*a*b_a, width=2*a, angle=theta, fill=False, color=color,lw=lw,ls='dotted')
    ax.add_artist(e)

masterlist = Table.read('table/AMUSING_PRAM_sample.v2.csv')
musecubes_path = "/net/ASTRO/pram/pram/starlight_el/bin1.v2.cA.CB16_6x16/"

fig, axs = plt.subplots(nrows = 8 , ncols = 6, figsize=(8,10))
axs = axs.ravel()
for i, name in enumerate(masterlist['Name']): #loops que passa por todos os arquivos 
    print(f"lendo a {name}")
    m = masterlist[i]
    cubegal = f'muse-{name}.bin1.v2.cA.CB16_6x16.dobby.fits'
    c=FitsCube(musecubes_path+cubegal) #carregando dados das galaxias

    x0, y0 = find_peak(x0=30, y0=30, image=c.fobs_norm)
    b_a = np.cos(np.radians(m['inc']))
    r50 = m['Reff']
    PA = m['PA']
    print(x0, y0, b_a, PA, r50)

    axs[i].imshow(c.fobs_norm)
    draw_r50(axs[i], x0, y0, r50, b_a, PA, color='w')
    
plt.show()  

   
