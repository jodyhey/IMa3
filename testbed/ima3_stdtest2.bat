ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8  -ha0.95 -hb0.9 -jhd
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -ha0.95 -hb0.9 -fima3_test_28.out -jhd
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_30.out  -b100000 -L5000  -j2 -c3 -g ima3_test_27.out.imapriors.txt
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_31.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_32.out -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_33.out -q10 -m1 -t1.5 -b100000 -L5000  -j03
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_34.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_35.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_34.out.mcf  -jh2 
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_36.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_37.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -ha0.97 -hb0.85 -jhd
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_38.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn24 -ha0.97 -hb0.85 -jhd
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_39.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_38.out.mcf   -hn24 -ha0.97 -jhd
mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -ha0.97 -c4
mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_41.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn20 -ha0.97 -c4 -fima3_test_40.out 
mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -ha0.97 -c4 -x 0 1 1e2 -x 2 3 0.1
