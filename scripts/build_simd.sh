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
gcc -c -g -Wall -O2 -fpic -o $FPIC_DIR/simd/scomplex.o $SRC_DIR/simd/complex_simd.c -lm
gcc -shared -o $LIBS_DIR/libscomplex.so fpic/simd/scomplex.o -lm

if [ $? -eq 0 ]; then
	echo "OK"

else
	echo "FAIL"
	rm -rf $FPIC_DIR
	exit $?

fi

rm -rf $FPIC_DIR
gcc -g -Wall -Llibs -Wl,-rpath=$LIBS_DIR $TEST_DIR/simd_test.c -lscomplex