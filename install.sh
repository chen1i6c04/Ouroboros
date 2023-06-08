#!/bin/bash -l

wget https://github.com/alekseyzimin/masurca/releases/download/v4.1.0/MaSuRCA-4.1.0.tar.gz
tar -xf MaSuRCA-4.1.0.tar.gz
rm MaSuRCA-4.1.0.tar.gz
mv MaSuRCA-4.1.0 MaSuRCA
cd MaSuRCA
./install.sh
