FPIC_DIR=fpic
LIBS_DIR=libs
SRC_DIR=src

# create the fpic directory
mkdir $FPIC_DIR

if [ ! -d "$LIBS_DIR" ];
then
	echo "Creating libs directory"
	mkdir $LIBS_DIR
fi

echo "Building and linking complex module"
gcc -c -g -fpic -o fpic/complex.o src/complex.c -lm
gcc -shared -o libs/libcomplex.so fpic/complex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	exit $?

fi

echo "Building and linking fft"
gcc -c -g -fpic -o $FPIC_DIR/fft.o $SRC_DIR/fft.c -lcomplex -lm
gcc -shared -o $LIBS_DIR/libfft.so $FPIC_DIR/fft.o $FPIC_DIR/complex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	exit $?

fi

echo "Building and linking recursive fft"
gcc -c -g -fpic -o $FPIC_DIR/rfft.o $SRC_DIR/recursive_fft.c -lfft -lcomplex -lm
gcc -shared -o $LIBS_DIR/librfft.so $FPIC_DIR/rfft.o $FPIC_DIR/fft.o $FPIC_DIR/complex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	exit $?

fi

echo "Building and linking iterative fft"
gcc -c -g -fpic -o $FPIC_DIR/itfft.o $SRC_DIR/iterative_fft.c -lfft -lcomplex -lm
gcc -shared -o $LIBS_DIR/libitfft.so $FPIC_DIR/itfft.o $FPIC_DIR/fft.o $FPIC_DIR/complex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	exit $?

fi
# echo "Building recursive fft"
# gcc -c -g -fpic -o $FPIC_DIR/ifft.o $SRC_DIR/iterative_fft.c -lfft -lcomplex -lm
# gcc -shared -o $LIBS_DIR/libifft.so $FPIC_DIR/ifft.o $FPIC_DIR/fft.o $FPIC_DIR/complex.o -lm

echo "Building and linking conv"
gcc -c -g -fpic -o $FPIC_DIR/conv.o $SRC_DIR/conv.c -litfft -lfft -lcomplex -lm
gcc -shared -o $LIBS_DIR/libconv.so $FPIC_DIR/conv.o $FPIC_DIR/itfft.o $FPIC_DIR/fft.o $FPIC_DIR/complex.o -lm

echo "Building and linking test"
gcc -g -Wall -Llibs -Wl,-rpath=$LIBS_DIR $SRC_DIR/test.c -lconv -litfft -lfft -lcomplex -lm

# remove the fpic directory
rm -r $FPIC_DIR

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	exit $?

fi

echo "Finished with 0 errors"
