#!/bin/bash

zcat $1 | grep -P 'taxon:(\d+)' -m 1 | grep -P '\d+' -o