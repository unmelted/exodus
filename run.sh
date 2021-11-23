rm libexodus.dylib
rm -rf generated/CMakeFiles
rm generated/cmake_install.cmake
rm generated/CMakeCache.txt

make -j4

FILE=libexodus.dylib
if [ -f "$FILE" ];then
    mv CMakeFiles generated/
    mv cmake_install.cmake generated/
    mv CMakeCache.txt generated/

    echo "Library created."
    cd streamlit
    streamlit run app.py
fi

if [ ! -d "generated" ];then
    mkdir generated
fi

if [ ! -d "log" ];then
    mkdir log
fi

if [ ! -d "saved" ];then
    mkdir log
fi