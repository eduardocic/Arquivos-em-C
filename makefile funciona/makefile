# =================================================================
# A referência empregada para escrever makefiles veio deste link:
# 
#   https://embarcados.com.br/introducao-ao-makefile/
# =================================================================
all: main

# Para criar o arquivo teste, eu preciso de?
main: teste.o main.o
	gcc -o main teste.o main.o

# Criando um objeto para 'main.c', o que eu preciso?
main.o: main.c teste.h
	gcc -o main.o main.c -c -W -Wall -ansi -pedantic
	
# Para eu criar um arquivo teste.o, o que eu preciso?	
teste.o: teste.c teste.h
	gcc -o teste.o teste.c -c -W -Wall -ansi -pedantic
	
clean:
	rm -rf *.o
