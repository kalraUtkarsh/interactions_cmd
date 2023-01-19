#!/usr/bin/env bash

python3.9 -m venv env

awk 'NF { print "export", $0}' configuration.env >> env/bin/activate

source env/bin/activate


# pip requirements
pip3 install wheel
pip3 install --upgrade setuptools
pip3 install -r pip-requirements.txt

