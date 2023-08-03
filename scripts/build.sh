SOURCE=${BASH_SOURCE[0]}
while [ -L "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )
  SOURCE=$(readlink "$SOURCE")
  [[ $SOURCE != /* ]] && SOURCE=$DIR/$SOURCE # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR=$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )

echo $DIR

FPIC_DIR="${DIR}/../fpic"
LIBS_DIR="${DIR}/../libs"
SRC_DIR="${DIR}/../src"

echo $LIBS_DIR

# create the fpic directory
mkdir $FPIC_DIR

if [ ! -d "$LIBS_DIR" ];
then
	echo "Creating shared directory"
	mkdir $LIBS_DIR
fi

echo "Building and linking complex module"
gcc -c -g -Wall -fpic -o $FPIC_DIR/complex.o src/complex.c -lm
gcc -shared -o $LIBS_DIR/libcomplex.so fpic/complex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -r $FPIC_DIR
	exit $?

fi

echo "Building and linking fft"
gcc -c -g -Wall -fpic -o $FPIC_DIR/fft.o $SRC_DIR/fft.c -lcomplex -lm
gcc -shared -o $LIBS_DIR/libfft.so $FPIC_DIR/fft.o $FPIC_DIR/complex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -r $FPIC_DIR
	exit $?

fi

# echo "Building and linking recursive fft"
# gcc -c -g -Wall -fpic -o $FPIC_DIR/rfft.o $SRC_DIR/recursive_fft.c -lfft
# gcc -shared -o $LIBS_DIR/librfft.so $FPIC_DIR/rfft.o $FPIC_DIR/fft.o

# if [ $? -eq 0 ]; then
# 	echo "OK"

# else
# 	echo "FAIL"
# 	rm -r $FPIC_DIR
# 	exit $?

# fi

echo "Building and linking iterative fft"
gcc -c -g -Wall -fpic -o $FPIC_DIR/itfft.o $SRC_DIR/iterative_fft.c -lfft
gcc -shared -o $LIBS_DIR/libitfft.so $FPIC_DIR/itfft.o $FPIC_DIR/fft.o

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -r $FPIC_DIR
	exit $?

fi
# echo "Building recursive fft"
# gcc -c -g -fpic -o $FPIC_DIR/ifft.o $SRC_DIR/iterative_fft.c -lfft -lcomplex -lm
# gcc -shared -o $LIBS_DIR/libifft.so $FPIC_DIR/ifft.o $FPIC_DIR/fft.o $FPIC_DIR/complex.o -lm

echo "Building and linking conv"
gcc -c -g -fpic -Wall -o $FPIC_DIR/conv.o $SRC_DIR/conv.c -litfft -lfft
gcc -shared -o $LIBS_DIR/libconv.so $FPIC_DIR/conv.o $FPIC_DIR/itfft.o $FPIC_DIR/fft.o

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -r $FPIC_DIR
	exit $?

fi

echo "Building and linking iterative fft radix 4"
gcc -c -g -Wall -fpic -o $FPIC_DIR/itfft4.o $SRC_DIR/fft_radix_4.c -lfft
gcc -shared -o $LIBS_DIR/libitfft4.so $FPIC_DIR/itfft4.o $FPIC_DIR/fft.o 

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -r $FPIC_DIR
	exit $?

fi

echo "Building and linking test"
gcc -g -Wall -Llibs -Wl,-rpath=$LIBS_DIR $SRC_DIR/test.c -lconv -litfft4 -litfft -lfft -lcomplex -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -r $FPIC_DIR
	exit $?

fi

# remove the fpic directory
rm -r $FPIC_DIR

# if [ $? -eq 0 ]; then
# 	echo "OK"

# else
# 	echo "FAIL"
# 	exit $?

# fi

echo "Finished with 0 errors"
