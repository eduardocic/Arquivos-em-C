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
tabela = pd.read_excel("lsim.xlsx")
#print(tabela)
tabela.columns= ['t', "y", "y"]

t   = tabela['t']
xk1 = tabela['xk1']
xk2 = tabela['xk2']
xh1 = tabela['xh1']
xh2 = tabela['xh2']
xh3 = tabela['xh3']
p11 = tabela['p11']
p22 = tabela['p22']
p33 = tabela['p33']

# plot do erro em altitude
err_x = xk1 - xh1
plt.plot(t, err_x, t, p11, '--k', t, -p11, '--k')
plt.xlabel("t");
plt.ylabel("Erro em Posição");
plt.axis([0, 30, -1500, 1500])
plt.grid()
plt.show()

# plot de erro em velocidade
err_v = xk2 - xh2
plt.figure()
plt.plot(t, err_v, t, p22, '--k', t, -p22, '--k')
plt.xlabel("t");
plt.ylabel("Erro em Velocidade");
plt.axis([0, 30, -500, 500])
plt.grid()
plt.show()


# plot de erro em aceleração
err_a = -32.2 - xh3;
plt.figure()
plt.plot(t, err_a, t, p33, '--k', t, -p33, '--k')
plt.xlabel("t");
plt.ylabel("Erro em Aceleração");
plt.axis([0, 30, -200, 200])
plt.grid()
plt.show()
