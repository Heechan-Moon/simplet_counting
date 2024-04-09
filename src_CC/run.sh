make clean
make main

ks="4 5 6"
ss="100 1000 10000 100000"
trial=5
versions="1 2" #1: original SC3, 2: memory-efficient SC3

#datas="coauth-DBLP coauth-MAG-Geology coauth-MAG-History
#        congress-bills contact-high-school contact-primary-school
#        DAWN email-Enron email-Eu NDC-classes NDC-substances
#        tags-ask-ubuntu tags-stack-overflow
#        threads-ask-ubuntu threads-math-sx threads-stack-overflow"
datas="toy"

echo "START"

for version in ${versions[@]}
do 
    for k in ${ks[@]}
    do
        for data in ${datas[@]}
        do
            echo "k="${k}
            for s in ${ss[@]}
            do
                echo "data="${data} "sample="${s} "version="${version}
                OMP_NUM_THREADS=6 ./main ${k} ${s} ${trial} ${data} ${version}
            done
        done
    done
done

echo "END"