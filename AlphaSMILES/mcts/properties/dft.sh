#!/usr/bin/env bash

export g09root=/usr/local/gaussian
source $g09root/g09/bsd/g09.profile

g09 $1
