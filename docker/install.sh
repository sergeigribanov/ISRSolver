export SOURCE_DIR=$HOME/source
export BUILD_DIR=$HOME/build
mkdir $SOURCE_DIR
mkdir $BUILD_DIR
mkdir $PKG_DIR
cd $SOURCE_DIR
git clone https://github.com/nlohmann/json
git clone https://github.com/GooFit/Minuit2
git clone --branch latest-stable https://github.com/root-project/root.git root_src
git clone https://github.com/sergeigribanov/ISRSolver
cd $BUILD_DIR
mkdir json
cd json
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/json  -DCMAKE_CXX_STANDARD=11 $SOURCE_DIR/json
make -j8
make install
mkdir $BUILD_DIR/minuit2
cd $BUILD_DIR/minuit2
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/minuit2  -DCMAKE_CXX_STANDARD=11 $SOURCE_DIR/Minuit2
make -j8
make install
mkdir $BUILD_DIR/root
cd $BUILD_DIR/root
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/root -DCMAKE_CXX_STANDARD=11 $SOURCE_DIR/root_src
make -j8
make install
source $PKG_DIR/root/bin/thisroot.sh
mkdir $BUILD_DIR/ISRSolver
cd $BUILD_DIR/ISRSolver
cmake -DCMAKE_INSTALL_PREFIX=$PKG_DIR/ISRSolver -Dnlohmann_json_DIR=$PKG_DIR/json/lib/cmake/nlohmann_json -DCMAKE_CXX_STANDARD=11 $SOURCE_DIR/ISRSolver
make -j8
make install
rm -rf $SOURCE_DIR
cd $HOME
export NOTEBOOKS_PREFIX=$HOME/workdir/notebooks
mkdir -p $NOTEBOOKS_PREFIX
mv yadisk.py $NOTEBOOKS_PREFIX/.
cd $NOTEBOOKS_PREFIX
python yadisk.py https://disk.yandex.com/d/XEes97fR2OsAqQ -p $NOTEBOOKS_PREFIX
tar -xf data.tar.gz
chmod -R a=rX $NOTEBOOKS_PREFIX/data
