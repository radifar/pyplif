#!/bin/bash

if  `cat $HOME/.bashrc | grep 'pyplif'`; then
    echo "PyPLIF is already installed"
    exit
else
    cp -r pyplif $HOME/.pyplif
    echo 'export pyplif="$HOME/.pyplif/pyplif.py"' >> $HOME/.bashrc
    echo "PyPLIF successfully installed"
fi
