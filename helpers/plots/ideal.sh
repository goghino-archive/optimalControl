#!/bin/bash

# Takes the following input:
# Id Nodes Time Ipopt_iter Schur_iter
# 1 1 63.9515 16 47
# 2 2 31.9894 16 47
# 3 4 16.0483 16 47
# 4 8 8.51151 16 47
# 5 16 8.51151 16 47
#
# and creates ideal scaling out of it, that is
# takes Time column and replaces the values
# by value of 1 node example divided by number of nodes.
# This results into following output
# 1 1 63.9515 16 47
# 2 2 31.9758 16 47
# 3 4 15.9879 16 47
# 4 8 7.99394 16 47
# 5 16 3.99697 16 47
awk 'NR<=1 { next } NR==2{a=$3} { print $1, $2,  a/$2, $4, $5 }' $1
