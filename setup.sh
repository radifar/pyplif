#!/bin/bash

if  `cat $HOME/.bashrc | grep 'pyplif'`; then
    echo "PyPLIF is already installed"
    exit
else
    if [ -e $HOME/.pyplif ]; then
        rm -r $HOME/.pyplif
    fi
    cp -r pyplif $HOME/.pyplif
    echo "alias pyplif='$HOME/.pyplif/pyplif.py'" >> $HOME/.bashrc
    echo "PyPLIF successfully installed"
fi
