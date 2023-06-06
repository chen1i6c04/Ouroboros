#!/bin/bash -l

if [[ ! -e $(which conda) ]];then
echo "Conda not found, please install Miniconda or Anaconda!"
exit 1
fi

wget https://github.com/alekseyzimin/masurca/releases/download/v4.1.0/MaSuRCA-4.1.0.tar.gz
tar -xf MaSuRCA-4.1.0.tar.gz
rm MaSuRCA-4.1.0.tar.gz
cd MaSuRCA-4.1.0
./install.sh
mv ./bin ..
cd ..
rm -R MaSuRCA-4.1.0
