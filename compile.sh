#! /bin/bash

exit_code=0

if [ $# -ne 0 ] && [ "$1" == "clean" ]; then
    message=" Cleaning up the project! " 
    color_code="01;34m"
    msg_length=${#message}

    term_width=$(tput cols)

    padding=$(( (term_width - msg_length) / 2))

    if (( (term_width - msg_length) % 2 != 0 )); then
      padding_right=$((padding + 1))
    else
      padding_right=$padding
    fi

    left_padding=$(printf "%*s" "$padding" "" | tr " " "=")
    right_padding=$(printf "%*s" "$padding_right" "" | tr " " "=")

    printf "\n\e[${color_code}${left_padding}${message}${right_padding}\e[0m\n"
    rm -rf run/* bin/* obj/*.o obj/*.d
    exit 0
  else
    term_width=$(tput cols)
    message=" Compiling Started "
    msg_length=${#message}

    padding=$(( (term_width - msg_length) / 2))
    if (( (term_width - msg_length) % 2 != 0 )); then
      padding_right=$((padding + 1))
    else
      padding_right=$padding
    fi
    left_padding=$(printf "%*s" "$padding" "" | tr " " "=")
    right_padding=$(printf "%*s" "$padding_right" "" | tr " " "=")

    printf "\n\e[01;34m${left_padding}${message}${right_padding}\e[0m\n"

    printf "\n"
    figlet -c -w $(tput cols) -f smslant Lipmann-Schwinger Solver
    printf "\n"

    printf "\e[01;34m%*s\e[0m\n" "$term_width" "" | tr " " "="
    printf "\n "

    mkdir build
    cd build
    cmake ..
    make -j $(getconf _NPROCESSORS_ONLN) 
    
    exit_code=$?


    if [[ $exit_code != 0 ]] ; then
      message=" Compile failed with errors! "
      color_code="01;31m"
    else 
      message=" Compiled successfully without errors! " 
      color_code="01;32m"
    fi
    
    msg_length=${#message}

    padding=$(( (term_width - msg_length) / 2))

    if (( (term_width - msg_length) % 2 != 0 )); then
      padding_right=$((padding + 1))
    else
      padding_right=$padding
    fi

    left_padding=$(printf "%*s" "$padding" "" | tr " " "=")
    right_padding=$(printf "%*s" "$padding_right" "" | tr " " "=")

    printf "\n\e[${color_code}${left_padding}${message}${right_padding}\e[0m\n"


    cd ..

    rm -rf build
fi

exit $exit_code
