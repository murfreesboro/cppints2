#
#   simple Makefile to build the module of types for debugging purpose
#
NAME       = test
CC		     = g++
CFLAGS     = -Wall -O3
INCLUDE    = -I../include 
LIB        = -lboost_filesystem -lboost_system
OBJC       = test.o ../src/shell.o ../src/basis.o  ../src/basisutil.o \
				 ../src/inttype.o ../src/integral.o  ../src/shellquartet.o ../src/derivinfor.o 

ALL: $(OBJC) 
	$(CC) -o  $(NAME) $(OBJC) $(LIB)

$(OBJC): %.o:%.cpp
	$(CC) $(INCLUDE) $(CFLAGS) -c  $< -o $@

.PHONY: clean
clean:
	@cd ../src; rm -rf *.o
	@rm -rf $(NAME) test.o

