FF=g++
PROG_NAME=exe
##ZIP_NAME=archive.zip
LFLAGS=
FFLAGS=-std=c++11
SRC=main.cpp fonction.cpp
##RES_DIR=res

OBJ=$(SRC:%.cpp)

all : $(PROG_NAME)

$(PROG_NAME) : $(SRC)
	$(FF) $(FFLAGS) -o $@ $^ 

clean :
	rm $(PROG_NAME)
	##rm -rf $(RES_DIR)
	##mkdir $(RES_DIR)	

##zip :
##	zip $(ZIP_NAME) $(SRC)