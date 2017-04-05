PROGRAM=PROG
BINDIR = bin
SRC = src
OBJ = obj
INC = include

CC     = g++
CLINK  = $(CC) 
C_LIB  = -lm
CFLAGS = -std=c++14 -Wall -O3 -I${INC} -fopenmp
CLINKFLAGS= -O3 -fopenmp


OBJS = ${OBJ}/main.o \
	${OBJ}/CPUGlobalOptimization.o \

${BINDIR}/${PROGRAM}: crdir ${OBJS}
	${CLINK} ${CLINKFLAGS} -o ${BINDIR}/${PROGRAM} ${OBJS} ${C_LIB}

crdir:
	mkdir -p obj
	mkdir -p bin
	
${OBJ}/main.o: ${SRC}/main.cpp ${INC}/interval.h ${INC}/CPUGlobalOptimization.h
	$(CC) $(CFLAGS) -c ${SRC}/main.cpp -o ${OBJ}/main.o
${OBJ}/CPUGlobalOptimization.o: ${SRC}/CPUGlobalOptimization.cpp ${INC}/CPUGlobalOptimization.h 
	$(CC) $(CFLAGS) -c ${SRC}/CPUGlobalOptimization.cpp -o ${OBJ}/CPUGlobalOptimization.o


clean:
	rm -rf ${OBJ}
	rm -rf bin
cleanall:
	rm -f ${OBJ}/*.o ${BINDIR}/${PROGRAM}
