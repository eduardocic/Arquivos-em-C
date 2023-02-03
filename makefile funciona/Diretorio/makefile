# Info diversas para facilitar o makefile
TARGET    = main
CC        = gcc
CFLAGS    = -c -W -Wall -ansi -pedantic -lm -std=c99
RM        = rm -rf

# Arquivos, headers e aqueles gerados do tipo '*.o'
 
C_SOURCE  = $(wildcard ./src/*.c)      
H_SOURCE  = $(wildcard ./include/*.h) 
O_FILES   = $(patsubst ./src/%.c, %.o, $(patsubst ./src,./include,$(C_SOURCE)))
TARGET_O  = $(TARGET).o
OBJ       = ./obj


# ======================================================= #
# 														  #	
# Executa o makefile, criando um arquivo do tipo 'main'   #
#  														  #
# ======================================================= #
all: obj main
	mv $(O_FILES)  $(OBJ)
	mv $(TARGET_O) $(OBJ)

main: $(O_FILES) $(TARGET_O)
	$(CC) -o $@ $(O_FILES) $(TARGET_O)

main.o: main.c $(H_SOURCE)
	$(CC) -I./include -o $@ main.c $(CFLAGS)

%.o: $(C_SOURCE) $(H_SOURCE)
	$(CC) -I./include -I./src -o $@ $(C_SOURCE) $(CFLAGS)

obj:
	mkdir $(OBJ)	
	
clean:
	$(RM) *.o
	$(RM) $(OBJ)/*.o
	rm $(TARGET)	
	$(RM) $(OBJ)
