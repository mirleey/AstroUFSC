from pycasso2 import FitsCube
import numpy as np
import matplotlib.pyplot as plt
from pycasso2.geometry import get_image_distance
from pycasso2.resampling import bin_edges
from astropy.table import Table
from pycasso2.geometry import find_peak



def draw_r50(ax, x0, y0, a, b_a, PA, color='indigo', lw=1.5):
    from matplotlib.patches import Ellipse
    center = np.array([x0, y0])
    #theta = PA - 90
    theta = PA
    e = Ellipse(center, height=2*a*b_a, width=2*a, angle=theta, fill=False, color=color,lw=lw,ls='dotted')
    ax.add_artist(e)


#plt.ion()
masterlist = Table.read('table/AMUSING_PRAM_sample.v2.csv')[:3]
musecubes_path = "/net/ASTRO/pram/pram/starlight_el/bin1.v2.cA.CB16_6x16/"
for i, name in enumerate(masterlist['Name']): #loops que passa por todos os arquivos 
    print(f"lendo a {name}")
    ml = masterlist[i]
    cubegal = f'muse-{name}.bin1.v2.cA.CB16_6x16.dobby.fits'
    c=FitsCube(musecubes_path+cubegal) #carregando dados das galaxias

    #m = (c.SN_normwin < 5) | (np.log10(c.EL_flux(6563)) < 2)
    x0, y0 = find_peak(x0=30, y0=30, image=c.fobs_norm, delta=15)
    b_a = np.cos(np.radians(ml['inc']))
    r50 = ml['Reff'] / c.pixelScale_arcsec
    r_max_pix = 4 * ml['Reff']
    PA = ml['PA']
    print(x0, y0, b_a, PA, r50)

    # FIXME: Melhorar essas máscara!
    m = (c.adev > 20) | (c.EL_flux(6563) <= 0)
    #m = (c.adev > 20) & (c.EL_flux(6563) < 1e-17)
    #m = (c.adev > 20) &(c.EL_flux(6563) <= 0)
    flux_ha = np.ma.array(np.log10(c.EL_flux(6563)), mask=m, copy=False)
    flux_n2 = np.ma.array(np.log10(c.EL_flux(6584)), mask=m, copy=False)
    n2ha = flux_n2 - flux_ha
    ewha = np.ma.array(c.EL_EW(6563), mask=m, copy=False)

    SF = ~m & (ewha > 3) & (n2ha < -0.4)
    AGN = ~m & (ewha > 6) & (n2ha >= -0.4)
    LINER = ~m & (ewha <= 6) & (ewha > 3) & (n2ha >= -0.4)
    HOLMES = ~m & (ewha <= 3)
    
    fig, axes = plt.subplots(nrows=3, ncols=2, layout= "constrained", figsize=(8, 9))
    fig.suptitle(name)

    ax = axes[0,0]
    im = ax.imshow(np.log10(c.fobs_norm), cmap="inferno_r")
    ax.scatter([x0],[y0], marker = "+", color = "r")
    fig.colorbar(im, ax=ax)
    font1 = {'family':'serif','color':'k','size':20}
    ax.set_title('Flux cont.', fontdict=font1)
    draw_r50(ax, x0, y0, r50, b_a, PA, color='limegreen')
    draw_r50(ax, x0, y0, 2 * r50, b_a, PA, color='limegreen')

    ax = axes[1,0]
    im = ax.imshow(np.log10(flux_ha), cmap="inferno_r")
    ax.scatter([x0],[y0], marker = "+", color = "r")
    fig.colorbar(im, ax=ax)
    font1 = {'family':'serif','color':'k','size':20}
    ax.set_title('Flux Ha', fontdict=font1)
    draw_r50(ax, x0, y0, r50, b_a, PA, color='r')
    draw_r50(ax, x0, y0, 2 * r50, b_a, PA, color='r')

    ax = axes[2,0]
    im = ax.imshow((ewha), cmap="inferno_r", vmax=6)
    ax.scatter([x0],[y0], marker = "+", color = "r")
    fig.colorbar(im, ax=ax)
    font1 = {'family':'serif','color':'k','size':20}
    ax.set_title('EW Ha', fontdict=font1)
    draw_r50(ax, x0, y0, r50, b_a, PA, color='r')
    draw_r50(ax, x0, y0, 2 * r50, b_a, PA, color='r')


#WHAN
    ax = axes[0,1]  
    ax.axvline(x=-0.4, color='k') # linha vertical
    ax.axhline(y=6, color='k')    # linha horizontal
    ax.axhline(y=3, color='k')
    ax.scatter(n2ha[SF], ewha[SF], marker='*', color='#ec407a')
    ax.scatter(n2ha[AGN], ewha[AGN], marker='x', color='steelblue')
    ax.scatter(n2ha[LINER], ewha[LINER], marker='D', color='#ffa726')
    ax.scatter(n2ha[HOLMES], ewha[HOLMES], marker='.', color='teal')
    font1 = {'family':'serif','color':'k','size':20}
    ax.set_title('WHaN', fontdict=font1)
    
#perfil radialN2Ha
    r = get_image_distance(flux_ha.shape, x0, y0, pa=np.radians(PA), ba=b_a)
    m = r < (r_max_pix)
    ax = axes[1,1]
    ax.axhline(y=-0.4, color='k') # linha vertical
    font1 = {'family':'serif','color':'k','size':20}
    ax.set_title('Perfil radial N2Ha', fontdict = font1)
    ax.scatter(r[m & SF], n2ha[m & SF], marker='*', color='#ec407a')
    ax.scatter(r[m & AGN], n2ha[m & AGN], marker='x', color='#7e57c2')
    ax.scatter(r[m & LINER], n2ha[m & LINER], marker='D', color='#ffa726')
    ax.scatter(r[m & HOLMES], n2ha[m & HOLMES], marker='.', color='teal')

    ax = axes[2,1]
    ax.scatter(r[m & SF], ewha[m & SF], marker='*', color='#ec407a')
    ax.scatter(r[m & AGN], ewha[m & AGN], marker='x', color='#7e57c2')
    ax.scatter(r[m & LINER], ewha[m & LINER], marker='D', color='#ffa726')
    ax.scatter(r[m & HOLMES], ewha[m & HOLMES], marker='.', color='teal')
    bins_r = np.arange(0, r_max_pix, 1)
    bins_e = bin_edges(bins_r)
    ew_r = c.radialProfile(ewha, bins_e, x0=x0, y0=y0, pa=np.radians(PA), ba=b_a, exact=False, mode="mean")
    ax.axhline(y=6, color='k')    # linha horizontal
    ax.axhline(y=3, color='k')
    #ax.axhline(y=1, color='indigo')
    #ax.axhline(y=2, color='indigo')
    ax.plot(bins_r, ew_r, color="indigo")
    font1 = {'family':'serif','color':'k','size':20}
    ax.set_title('Perfil Radial EW Ha', fontdict = font1)
    ax.set_ylim(0,25)
    plt.savefig(f'radprof/{name}.png')



#plt.show()
