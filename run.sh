#!/bin/bash
make clean 
make -j8
for i in `seq 1 10`;
 do
     ./main 0405e mc_prph 0 0 0 $i
 done   


#for i in `seq 1 10`;
# do
#     $i
# done
rm mc_prph0405e_all.root
hadd mc_prph0405e_all.root mc_prph0405e_parton_10.root mc_prph0405e_parton_9.root mc_prph0405e_parton_8.root mc_prph0405e_parton_7.root mc_prph0405e_parton_6.root mc_prph0405e_parton_5.root mc_prph0405e_parton_4.root mc_prph0405e_parton_3.root mc_prph0405e_parton_2.root mc_prph0405e_parton_1.root

echo "your computation is done" | mail -s "ready" "greyxray@gmail.com"
#hadd mc_prph0405e_all.root mc_prph0405e_parton10.root mc_prph0405e_parton9.root mc_prph0405e_parton8.root mc_prph0405e_parton7.root mc_prph0405e_parton6.root mc_prph0405e_parton5.root mc_prph0405e_parton4.root mc_prph0405e_parton3.root mc_prph0405e_parton2.root mc_prph0405e_parton1.root
