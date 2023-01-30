import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


tabela = pd.read_excel("random.xlsx")
#print(tabela)

plt.plot(tabela['t'], tabela['x_real'], tabela['t'], tabela['x_com_ruido'], 'r')
plt.xlabel("t")
plt.ylabel("x")
plt.grid()
plt.show()



plt.figure();
plt.plot(tabela['t'], tabela['v_real'])
plt.xlabel("t")
plt.ylabel("v")
plt.grid()
plt.show()
