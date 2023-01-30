import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# Para executar o arquivo em C primeiro.
#import subprocess        
#
# Arquivos em C os quais eu preciso rodar para finalizar a execução      
#subprocess.call(["gcc", "main2.c", "bib_matriz.c"]) 
#tmp = subprocess.call("./main")


#if temp = 0
#tabela = pd.read_table("test.txt", delimiter='\n')
tabela = pd.read_excel("hist.xlsx")
print(tabela)
tabela.columns= ['t']

t   = tabela['t']

plt.histogram(t)
plt.grid()
plt.show()

