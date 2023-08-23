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
TEST_DIR="${DIR}/../test"

# create the fpic directory
mkdir $FPIC_DIR
mkdir $FPIC_DIR/simd

if [ ! -d "$LIBS_DIR" ];
then
	echo "Creating shared directory"
	mkdir $LIBS_DIR
fi

echo "Building and linking complex simd module"
gcc -c -g -Wall -fpic -o $FPIC_DIR/simd/scomplex.o $SRC_DIR/simd/complex_simd.c -lm
gcc -shared -o $LIBS_DIR/libscomplex.so $FPIC_DIR/simd/scomplex.o -lm

echo "Building and linking fft simd module"
gcc -c -g -Wall -fpic -o $FPIC_DIR/simd/sfft.o $SRC_DIR/simd/fft_simd.c -lscomplex -lm
gcc -shared -o $LIBS_DIR/libsfft.so $FPIC_DIR/simd/sfft.o $FPIC_DIR/simd/scomplex.o -lm

echo "Building and linking fftr2 module"
gcc -c -g -Wall -fpic -o $FPIC_DIR/simd/sfftr2.o $SRC_DIR/simd/simd_fft_2.c -lsfft -lscomplex -lm
gcc -shared -o $LIBS_DIR/libsfftr2.so $FPIC_DIR/simd/sfftr2.o $FPIC_DIR/simd/sfft.o $FPIC_DIR/simd/scomplex.o -lm

echo "Building and linking fftr4 module"
gcc -c -g -Wall -fpic -o $FPIC_DIR/simd/sfftbr4.o $SRC_DIR/simd/simd_fft_br4.c -lsfft -lscomplex -lm
gcc -shared -o $LIBS_DIR/libsfftbr4.so $FPIC_DIR/simd/sfftbr4.o $FPIC_DIR/simd/sfft.o $FPIC_DIR/simd/scomplex.o -lm

gcc -c -g -Wall -fpic -o $FPIC_DIR/simd/sfftr4.o $SRC_DIR/simd/simd_fft_4.c -lsfftbr4 -lsfft -lscomplex -lm
gcc -shared -o $LIBS_DIR/libsfftr4.so $FPIC_DIR/simd/sfftr4.o $FPIC_DIR/simd/sfftbr4.o $FPIC_DIR/simd/sfft.o $FPIC_DIR/simd/scomplex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -rf $FPIC_DIR
	exit $?

fi

rm -rf $FPIC_DIR
gcc -g -Wall -Llibs -Wl,-rpath=$LIBS_DIR $TEST_DIR/simd_test.c -o simd_test.out -lscomplex
gcc -g -Wall -Llibs -Wl,-rpath=$LIBS_DIR $TEST_DIR/fft_twiddle_test.c -o fft_twid_test.out -lscomplex -lsfft
gcc -g -Wall -Llibs -Wl,-rpath=$LIBS_DIR $TEST_DIR/fft_r2_test.c -o fft_r2_test.out -lsfftr2 -lsfft -lscomplex 
gcc -g -Wall -Llibs -Wl,-rpath=$LIBS_DIR $TEST_DIR/fft_r4_test.c -o fft_r4_test.out -lsfftr4 -lsfftbr4 -lsfft -lscomplex 