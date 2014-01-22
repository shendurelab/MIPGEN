OBJS = mipgen.o Featurev5.o SVMipv4.o PlusSVMipv4.o MinusSVMipv4.o svm.o
CC = g++
DEBUG = -g
CFLAGS = -Wall -c ${DEBUG}
LFLAGS = -Wall ${DEBUG}

mipgen : ${OBJS}
	${CC} ${LFLAGS} ${OBJS} -o mipgen  
mipgen.o : mipgen.cpp
	${CC} ${CFLAGS} mipgen.cpp -I ./
SVMipv4.o : SVMipv4.cpp SVMipv4.h
	${CC} ${CFLAGS} SVMipv4.cpp -o SVMipv4.o
PlusSVMipv4.o : PlusSVMipv4.cpp PlusSVMipv4.h
	${CC} ${CFLAGS} PlusSVMipv4.cpp -o PlusSVMipv4.o 
MinusSVMipv4.o : MinusSVMipv4.cpp MinusSVMipv4.h
	${CC} ${CFLAGS} MinusSVMipv4.cpp -o MinusSVMipv4.o
svm.o : svm.cpp svm.h
	${CC} ${CFLAGS} svm.cpp -o svm.o
Featurev5.o : Featurev5.cpp Featurev5.h
	${CC} ${CFLAGS} Featurev5.cpp -o Featurev5.o
clean :
	\rm *.o
