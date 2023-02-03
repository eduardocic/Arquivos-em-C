# Info diversas para facilitar o makefile
TARGET    = main
TARGET_C  = $(TARGET).c
CC        = gcc
CFLAGS    = -lm -std=c99 #-W -Wall -ansi -pedantic 
RM        = rm -rf

# Arquivos, headers e aqueles gerados do tipo '*.o'
 
C_SOURCE  = $(wildcard ./src/*.c)      
H_SOURCE  = $(wildcard ./include/*.h) 
DOT_C     = $(patsubst ./src/%.c, %.c, $(patsubst ./src,./include,$(C_SOURCE)))
DOT_H     = $(patsubst ./include/%.h, %.c, $(patsubst ./src,./include,$(H_SOURCE)))
OBJ       = ./obj


# ======================================================= #
# 														  #	
# Executa o makefile, criando um arquivo do tipo 'main'   #
#  														  #
# ======================================================= #
all:
	mv ./src/*.c .
	mv ./include/*.h .	
	@ echo $(DOT_C)
	@ echo $(DOT_H)
	$(CC) -o $(TARGET) $(TARGET_C) $(DOT_C) $(CFLAGS)
	mv *.c ./src/
	mv *.h ./include/
	mv ./src/$(TARGET_C) .

clean:
	$(RM) $(TARGET)
	mv *.c ./src/
	mv *.h ./include/
	mv ./src/$(TARGET_C) .
	
.PHONY: all clean	
