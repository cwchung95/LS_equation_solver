#! /bin/bash

exit_code = 0;

if [ $# -ne 0 ] && [ "$1" == "clean" ]; then
    echo "Cleaning up..."
    rm -rf run/* bin/* obj/*.o obj/*.d
    exit 0
  else
    mkdir build
    cd build
    cmake ..
    make -j $(getconf _NPROCESSORS_ONLN) 
    exit_code = $?

    if [[ $exit_code != 0 ]] ; then
        echo "\n\e[01;31m===== ERRORS and WARNINGS =====\e[0m\n"
        exit $exit_code
    else 
      printf "\n\e[01;32m Compiled successfully without errors or warnings! \e[0m\n\n"
    fi

    cd ..

    rm -rf build
fi

exit $exit_code
