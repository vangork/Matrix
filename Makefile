GCC = g++
MAT = Matrix
LIBNAME = libmatly
TARGET = test

${TARGET}: ./src/${TARGET}.cpp ${LIBNAME}.a
	${GCC} -c -pipe -O2 -fPIC -I./include -o ./src/${TARGET}.o ./src/${TARGET}.cpp
	${GCC} -o ${TARGET} ./lib/${LIBNAME}.a ./src/${TARGET}.o -llapack   

${LIBNAME}.a: ${MAT}.o
	rm -f ./lib/${LIBNAME}.a
	ar cqs ${LIBNAME}.a ./src/${MAT}.o
	mv ${LIBNAME}.a ./lib/

${MAT}.o: ./src/${MAT}.cpp
	${GCC} -c -pipe -O2 -fPIC -I./include -o ./src/${MAT}.o ./src/${MAT}.cpp

clean:
	rm ./src/*.o ${TARGET} ./lib/${LIBNAME}.a -f

