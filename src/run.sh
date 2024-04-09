make clean
make main

ks="4 5 6"
ss="100 1000 10000 100000"
try="5"
threads=6

#datas="coauth-DBLP coauth-MAG-Geology coauth-MAG-History
#        congress-bills contact-high-school contact-primary-school
#        DAWN email-Enron email-Eu NDC-classes NDC-substances
#        tags-ask-ubuntu tags-stack-overflow
#        threads-ask-ubuntu threads-math-sx threads-stack-overflow"
datas="toy"

echo "START"

for k in ${ks[@]}
do
    for data in ${datas[@]}
    do
        for s in ${ss[@]}
        do
            echo "data="${data} "k="${k} "sample="${s} "try="${try} 
            OMP_NUM_THREADS=${threads} ./main ${data} ${k} ${s} ${try} 
        done
    done
done

echo "END"
