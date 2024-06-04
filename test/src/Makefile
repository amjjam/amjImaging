FLAGS=-g -Wall -I /usr/lib/x86_64-linux-gnu/openmpi/include \
    -I ../../include -I /opt/amj/include -pg -std=c++11 
CPPFLAGS=$(FLAGS)
CPP=g++
MPICPP=mpic++

TEST1_SRC=test1.C
TEST1_OBJ=$(TEST1_SRC:.C=.o)

OBJS=$(TEST1_OBJ)

BINS=test1


all: $(BINS)

test1: $(TEST1_OBJ)
	$(MPICPP) -o $@ $^ -L/opt/amj/lib -laVec -lm -lrt -lmpi -lgsl
clean:
	- rm $(BINS)
	- rm $(OBJS)
