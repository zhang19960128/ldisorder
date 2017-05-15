#!/bin/bash
icpc -mkl=sequential -std=c++11 -o go function.cpp main.cpp;
for ita in 0.0 0.6 1.0 1.4 1.56 1.72 1.80 1.88 1.92 1.96
do
cat>run${ita}.sh<<EOF
#!/bin/bash
./go ${ita} eigen${ita}.txt;
EOF
chmod +x run${ita}.sh
bsub -q "short" ./run${ita}.sh
done
