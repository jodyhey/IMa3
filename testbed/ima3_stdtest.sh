./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4	
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -hfg -ha0.95 -hb0.8 -d10 -r1 -c3  -gima3_priorfile_3pops.txt
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_7.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hfs -ha0.98 -d5
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_8.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j5 -p4
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_9.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j6 -p4
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_10.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j1 -p4
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_11.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j7 -p4
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_12.out -q10 -t1.5 -b100000 -L5000  -r1 -j8 -p4
./IMa3_stdtest -help  > ima3_test_13.out 
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_14.out -q10 -m1 -t1.5 -b100000 -L5000  -r1  -j9
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_15.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -jx
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_16.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_17.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -c0 
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_18.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 
./IMa3_stdtest -i ima3test4pop5loci.u  -o ima3_test_19.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -r5
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_20.out  -q10 -m1 -t1.5 -b100000 -L5000  -j01
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -ha0.95 -hb0.9
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8 -ha0.95 -hb0.9
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -ha0.95 -hb0.9 -fima3_test_26.out
./IMa3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -b100000 -L5000  -j2 -c3 -g ima3_test_25.out.imapriors.txt
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_1.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
./IMa3_stdtest -i ima3test8pop2loci.u  -o ima3_test_30.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
./IMa3_stdtest -i ima3test8pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j03
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_jh2_3.out.mcf  -jh2 
./IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_16.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2
mpirun -n 4 IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_mpi_1.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -ha0.97 -hb0.85
mpirun -n 4 IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_mpi_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn12 -ha0.97 -hb0.85
mpirun -n 4 IMa3_stdtest -i ima3test3pop4loci.u -o ima3_test_mpi_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_mpi_3.out.mcf   -hn12 -ha0.97 -hb0.85
mpirun -n 4 IMa3_stdtest -i ima3test4pop2loci.u -o ima3_test_mpi_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8 -ha0.95 -hb0.9
mpirun -n 4 IMa3_stdtest -i ima3test4pop2loci.u -o ima3_test_mpi_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -ha0.95 -hb0.9 -fima3_test_mpi_26.out
