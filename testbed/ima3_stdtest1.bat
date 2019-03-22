ima3_stdtest -help  > ima3_test_0.out 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -ha0.95 -hb0.8 -d10 -r1 -c3  -gima3_priorfile_3pops.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_7.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -ha0.98 -d5 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_8.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -d5
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_9.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j5 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_10.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j6 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_11.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j1 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_12.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j7 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_13.out -q10 -t1.5 -b100000 -L5000  -r1 -j8 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_14.out -q10 -m1 -t1.5 -b100000 -L5000  -r1  -j9
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_15.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -jx
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_16.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_17.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -c0 
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_18.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 
ima3_stdtest -i ima3test4pop5loci.u  -o ima3_test_19.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -r5
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_20.out  -q10 -m1 -t1.5 -b100000 -L5000  -j01
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -ha0.95 -hb0.9
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -ha0.97 
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -c4
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123 -x 0 1 1e2 -x 2 3 0.1
