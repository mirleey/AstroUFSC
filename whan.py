from astropy.table import Table
import os
import matplotlib.pyplot as plt
plt.ion()           #MADIÇÃO interativo
from pathlib import Path
import numpy as np
from pycasso2 import FitsCube

cubenames=os.listdir("/net/ASTRO/manga/data/PRAM/muse/starlight_el/bin21.cA.CB16_6x16") #juntando todas os arquivos, se usar o (!ls) ativamos o linux, portanto não se usa (), "" etc.      
# all_galaxy=Path('/net/ASTRO/manga/data/PRAM/muse/starlight_el/bin21.cA.CB16_6x16/') #juntando dois arquivos
musecubes_path='/net/ASTRO/manga/data/PRAM/muse/starlight_el/bin21.cA.CB16_6x16/'
#galaxy=all_galaxy/'muse-MCG-01-38-015.bin21.cA.CB16_6x16.dobby.fits'             
# galaxy_one=all_galaxy/'muse-3C264.bin21.cA.CB16_6x16.dobby.fits'


sir_mixe= []
mixe_bq=[]

for cubegal in cubenames: #loops que passa por todos os arquivos
    
    #print(galaxy)


    c=FitsCube(musecubes_path+cubegal) #carregando dados das galaxias
    #c.f_obs
    #c.plotPixel(1,1) #desenho dos espectros dos modelos
    #plt.show() #medida de 1pixel

    #(log)ha é 6563 e (log)nII é 6584 #Fluxo
    fluxonii= c.EL_integ(6584)['EL_F']
    fluxoha=c.EL_integ(6563)['EL_F']

    #np.log10(fluxonii/fluxoha)#calcular o lognii/ha 
    n2ha=np.log10(fluxonii/fluxoha)
    larguraha=c.EL_integ(6563)['EL_EW'] #largura equivalente


    ewha=np.log10(larguraha) #posição mo gráfico no eixo de y

    sir_mixe.append(n2ha)
    mixe_bq.append(ewha) #ad os valores dentro de listas

    #plt.imshow(c.EL_flux(6563)) #imagem com fluxo de ha
    #c.EL_integ(6563)['EL_F'] # fluxo integrando a galaxia inteira em ha
    #c.EL_integ(6563)['EL_EW']#largura equivalente de ha
    #dtype=(numpy.record, [('lambda', '>i4'), ('line', 'S20'), ('El_l0', '>f8'), ('El_F', '>f8'), ('El_v0', '>f8'), ('El_vd', '>f8'), ('El_flag', '>i4'), ('El_EW', '>f8'), ('El_lcrms', '>f8'), ('El_vdins', '>f8')]))


plt.scatter(sir_mixe, mixe_bq) #plot dos pontos

#após plotar, vamos salvar os valores em tabelas
lapis_dasilva= Table() #criando tabela e guardando os valores em colunas
lapis_dasilva.add_column(sir_mixe, name='n2ha')
lapis_dasilva.add_column(mixe_bq, name='ewha')

#criando um txt e formatando dados da tabela e guardando as infos
lapis_dasilva.write('whan.txt', format='ascii.commented_header', overwrite=True)

