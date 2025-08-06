i=$1

cp coby_input.top coby_input.yaml  protonations.txt $i/

cp -r toppar/ Scripts/ $i/

cd $i/

python Scripts/AABY.py -f $i.pdb -r 4 --chains A,B,C,D --ssbond --protlist protonations.txt --coby-yaml coby_input.yaml

for rep in {1..4}; do
    cd r$rep
    for j in $(ls toppar/system*itp); do python ../../Scripts/generate_posres_full.py --itp $j; done
    j=CHL; python ../../Scripts/generate_posres_full.py --itp toppar/$j.itp

    python ../../Scripts/prepare_for_simulations.py ${i}_AABY
    cd ../
done

cd ../
#
#
#
#i=$1
#
#cp coby_input.top coby_input.yaml  protonations.txt $i/
#
#cp -r toppar/ Scripts/ $i/
#
#cd $i/
#
#python Scripts/AABY.py -f $i.pdb -r 4 --chains A,B,C,D --ssbond --protlist protonations.txt --coby-yaml coby_input.yaml
#
#for j in $(ls toppar/system*itp); do python ../Scripts/generate_posres_full.py --itp $j; done
#j=CHL; python ../Scripts/generate_posres_full.py --itp toppar/$j.itp
#
#python Scripts/prepare_for_simulations.py ${i}_AABY
#
#cd ../
