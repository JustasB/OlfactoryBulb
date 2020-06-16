#!/bin/bash

# Fail on first error
set -e

# Remove existing .mod file compilations
find . -name x86_64 -exec sudo rm -rif {} || true \;

# Recompile mod files of all previous models
cd BhallaBower1993/; nrnivmodl; cd ..;
cd Birgiolas2020/Mechanisms/; nrnivmodl; cd ../..;
cd Chen2002/; nrnivmodl; cd ..;
cd David2008/; nrnivmodl; cd ..;
cd Davison2000/; nrnivmodl; cd ..;
cd Davison2003/; nrnivmodl; cd ..;
cd Djurisic2008/; nrnivmodl; cd ..;
cd KaplanLansner2014/; nrnivmodl; cd ..;
cd LiCleland2013/; nrnivmodl; cd ..;
cd McTavish2012/; nrnivmodl; cd ..;
cd Migliore2005GJs/; nrnivmodl; cd ..;
cd Migliore2007columns/; nrnivmodl; cd ..;
cd Migliore2008microcircuits/; nrnivmodl; cd ..;
cd Migliore2014bulb3d/; nrnivmodl; cd ..;
cd Migliore2015operators/; nrnivmodl; cd ..;
cd MiglioreMcTavish2013/; nrnivmodl; cd ..;
cd Oconnor2012/; nrnivmodl; cd ..;
cd Popovic2005/; nrnivmodl; cd ..;
cd RubinCleland2006/; nrnivmodl; cd ..;
cd Saghatelyan2005/; nrnivmodl; cd ..;
cd Shen1999/; nrnivmodl; cd ..;
cd Short2016/; nrnivmodl; cd ..;
cd Yu2012/; nrnivmodl; cd ..;
