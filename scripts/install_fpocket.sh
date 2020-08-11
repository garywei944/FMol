#!/usr/bin/env bash
echo "
-------------------------------
----- Installing fpocket ------
-------------------------------
"

# BUG: may cannot access /tmp
cd /tmp || exit
sudo apt -y install libnetcdf-dev
git clone https://github.com/Discngine/fpocket.git
cd fpocket || exit
make -j "$(nproc)"
sudo make install
cd .. || exit
rm -fr fpocket

echo "
-------------------------------------------
----- Successfully installed fpocket ------
-------------------------------------------
"
