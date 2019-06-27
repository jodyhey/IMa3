NOTE:  the most current version of this and all testbed files are in E:\genemod\ML-MCMC\SEAI\IMa3\testbed
DO NOT use older versions of testbed files from old directories. 
THIS IS also true for the linux version IMa3_stdtest.sh
basic instructions:
	all the main files are in E:\genemod\ML-MCMC\SEAI\IMa3\testbed

	to do a run:
		copy the files from E:\genemod\ML-MCMC\SEAI\IMa3\testbed into a new directory
		e.g. E:\genemod\ML-MCMC\SEAI\IMa3_work\stdTest\####
		
		delete any previous results file in that folder
		
		compile a new IMa3_stdtest.exe file with STDTEST on 
		copy IMa3_stdtest.exe to the new directory 
		run ima3_stdtest1.bat  and ima3_stdtest2.bat
		
	for wsl
		from unbuntu
		copy the files from  E:\genemod\ML-MCMC\SEAI\IMa3\testbed into a new directory 
			e.g. E:\genemod\ML-MCMC\SEAI\IMa3_work\stdTest\wsl\####
			delete any previous results file in that folder
		compile a new IMa3_stdtest  using make testbed
		run ima3_stdtest.sh

	current standard is 3_5_2019

6/26/2019
	made some changes to prevent crash due to filled migration array
	stdtest results looked unchanged
	
	
3/25/2019
	the 3/21/2019 run looked good.  only test22 looked different from the 3/8/2019 runs
	
	design a new test22 that uses the sigmoid model with lots more chains 
	
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b30000  -L5000  -j013  -hn100  -ha0.999

	now set the swap interval to 5 steps 
		#define CHAINSWAPINTERVAL 4 
		
	and rerun the full set to set the new standard 
	
	ima3_stdtest1.bat:
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
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b30000  -L5000  -j013  -hn100  -ha0.999
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -c4
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123 -x 0 1 1e2 -x 2 3 0.1


	ima3_stdtest2.bat
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8  -ha0.95 -hb0.9
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -ha0.95 -hb0.9 -fima3_test_28.out
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_30.out  -b100000 -L5000  -j2 -c3 -g ima3_test_27.out.imapriors.txt
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_31.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
	ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_32.out -q10 -m1 -t1.5 -b100000 -L5000  -j013
	ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_33.out -q10 -m1 -t1.5 -b100000 -L5000  -j03
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_34.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_35.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_34.out.mcf  -jh2 
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_36.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2
	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_37.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -ha0.97 -hb0.85
	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_38.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn24 -ha0.97 -hb0.85
	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_39.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_38.out.mcf   -hn24 -ha0.97
	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -ha0.97 -c4
	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_41.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn20 -ha0.97 -c4 -fima3_test_40.out 
	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -ha0.97 -c4 -x 0 1 1e2 -x 2 3 0.1
	
	
	
		
3/21/2019

I changed the default mode of operation for multiple chains so that the prior is raised to beta by default
now you have to use -jhD to turn that off
also -c4  will turn that off 
this will affect lots of results with multiple chains and when -c4 is not used  
change commands for  6, 21, 28, 29, 37, 38, 39 to remove -jhd from the command line 

This will also cause test 22  to be different,  swapping rates should be much less and mixing should be poor

I also changed the interval between swap updates to 4 from 0 
	this will cause lots of differences
	for this run,  reset the interval to 0,  and check (to see effect of changing -jhD)
	then do another run with the new swap interval to set the new standard 

	ima3_stdtest1.bat:
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


	ima3_stdtest2.bat
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8  -ha0.95 -hb0.9
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -ha0.95 -hb0.9 -fima3_test_28.out
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_30.out  -b100000 -L5000  -j2 -c3 -g ima3_test_27.out.imapriors.txt
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_31.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
	ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_32.out -q10 -m1 -t1.5 -b100000 -L5000  -j013
	ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_33.out -q10 -m1 -t1.5 -b100000 -L5000  -j03
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_34.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_35.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_34.out.mcf  -jh2 
	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_36.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2
	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_37.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -ha0.97 -hb0.85
	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_38.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn24 -ha0.97 -hb0.85
	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_39.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_38.out.mcf   -hn24 -ha0.97
	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -ha0.97 -c4
	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_41.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn20 -ha0.97 -c4 -fima3_test_40.out 
	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -ha0.97 -c4 -x 0 1 1e2 -x 2 3 0.1



3/8/2019
	some smallish changes,  looks the same as 3/5/2019
	
3/5/2019
	did some more fixes to the code
		now the only jobs that turn up as different from 2/24/2019 set are 22, 23, 40 and 41 
		these are all -j013 jobs using no  or using -c4 
		
		jobs 7 and 8  which used no -c4 were fine 
		
		so the main problem seems to have been when using NOjhD and -j0 
		
	current list of jobs
	ima3_stdtest1.bat:
	ima3_stdtest -help  > ima3_test_0.out 
	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4
	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -ha0.95 -hb0.8 -jhd -d10 -r1 -c3  -gima3_priorfile_3pops.txt
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
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -ha0.95 -hb0.9 -jhd
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -ha0.97 
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -c4
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123 -x 0 1 1e2 -x 2 3 0.1


	ima3_stdtest2.bat
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
		
3/3/2019

found a significant bug that distorted posteriors of topologies.  very bad

on 5/9/2018 I switched the default heating scheme to not raising the prior to beta 
	i.e. lilihood^beta + prior,   rather than  (likelihood+prior)^beta
	this is the heating scheme required of thermodyanmic integration 
	
	at that point I set -jhD  to turn on the traditional heating scheme i.e. (likelihood+prior)^beta
	
	but I noticed that it was giving me weird tree for mus  
	
	
	found the likely cause in swap_weight()
		this should affect runs using multiple chains and that did not use -jhD 
		runs marked -jhD should be the same 

		which runs are multiple chains do not raise prior to beta? 

		22  is -hn12 and -j013 and no -jhD

		23 is -hn12 and -j013 and -c4  (no -jhD)

		40 is -j013,  -hn20  and -c4 

		41 is -j013 -hn20  and -c4  
		
		7, 8  used -c4 

	run the code with the bug fix in 3_3_2019
	
	affected jobs
	22, 41, 23, 40, 7, 8, 


1/22/2019

worked on handling small errors in input file format better 
also modified ima3test3pop4loci.u  to have various small formatting errors in gene names 

reran testbed compare with 8/20/2018
also run testbed on wsl linux 
looks ok



7/5/2018

added to main output file a table of estimate locations near the top 

reran testbed, looked ok

6/5/2018

made some changes:
	stop printing mutation scalar trend plots
	changed the way -x works and when it is used
		user can now specify any priors on clades
		also start with a sample of trees from the prior

add command line to stdtest1.bat 
	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123 -x 0 1 1e2 -x 2 3 0.1

add command line to stdtest2.bat 
	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40x.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -ha0.97 -c4 -x 0 1 1e2 -x 2 3 0.1

update 	ima3_stdtest.sh

run in E:\genemod\ML-MCMC\SEAI\IMa3_work\stdTest\6_5_2018



5/18/2018

made an ubuntu 18.04 (WSL) version 
	runs in  ~/IMa3/run/stdtest/
	uses   ima3_stdtest.sh
	
	ran it and copied all files to 
		E:\genemod\ML-MCMC\SEAI\IMa3_work\stdTest\ubuntu_5_18_2018
	so can compare ubuntu runs in the future 
	
	
5/16/2018

changed treestring format to do without colons

this means changing the way prior files are written and read 

so make a new version of ima3_priorfile_3pops.txt  without colons after parentheses

compare with 5/9/2018 set 

looks good 5/16/2018 is new standard 

copied all the standard files into 
E:\genemod\ML-MCMC\SEAI\IMa3\testbed



5/9/2018
got rid of -hm  (fomerly -hfs)
	now -ha only invokes HFULL
		-ha and -hb invokes HGEOMETRIC
		neither -ha nor -hb invokes HEVEN
	

updated contents of  ima3_stdtest1.bat:
ima3_stdtest -help  > ima3_test_0.out 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -ha0.95 -hb0.8 -jhd -d10 -r1 -c3  -gima3_priorfile_3pops.txt
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
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -ha0.95 -hb0.9 -jhd
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -ha0.97 
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -c4
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123

updated contents of ima3_stdtest2.bat
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
		
For all this editing
	5/9/2018  is the set to compare stdtest runs with 

5/8/2018
merged with Jared's version 
	brings in xml coding 
	first stdtest run with xml was on 5/8/2018
	compared 5/8/2018 run with 4/30/2018
		all matched (except timings)  and 
		still have a bug where autoc_vals[0] are getting overwritten 
			the first value of many of these is getting overwritten 
			
to compile for stdtest,  turn on 
	XMLOUTPUT
	STDTEST
	
two batch files:
ima3_stdtest1.bat  and  ima3_stdtest2.bat

contents of  ima3_stdtest1.bat:
ima3_stdtest -help  > ima3_test_0.out 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -hmg -ha0.95 -hb0.8 -jhd -d10 -r1 -c3  -gima3_priorfile_3pops.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_7.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hms -ha0.98 -d5 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_8.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hme -d5
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
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -hmg -ha0.95 -hb0.9 -jhd
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -hms -ha0.97 
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -hme -c4
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8 -hmg -ha0.95 -hb0.9 -jhd
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -hmg -ha0.95 -hb0.9 -fima3_test_28.out -jhd

contents of ima3_stdtest2.bat
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_30.out  -b100000 -L5000  -j2 -c3 -g ima3_test_27.out.imapriors.txt
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_31.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_32.out -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_33.out -q10 -m1 -t1.5 -b100000 -L5000  -j03
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_34.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_35.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_34.out.mcf  -jh2 
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_36.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_37.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -hmg -ha0.97 -hb0.85 -jhd
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_38.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn24 -hmg -ha0.97 -hb0.85 -jhd
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_39.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_38.out.mcf   -hn24 -hms -ha0.97 -jhd
mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -hms -ha0.97 -c4
mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_41.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn20 -hms -ha0.97 -c4 -fima3_test_40.out 




	
	

5/1/2018
old#   new command
13 	ima3_stdtest -help  > ima3_test_0.out 
1	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4
2	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
3	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
4	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
5	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
6	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -hmg -ha0.95 -hb0.8 -jhd -d10 -r1 -c3  -gima3_priorfile_3pops.txt
7	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_7.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hms -ha0.98 -d5 
7	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_8.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hme -d5
8	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_9.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j5 -p4
9	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_10.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j6 -p4
10	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_11.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j1 -p4
11	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_12.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j7 -p4
12	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_13.out -q10 -t1.5 -b100000 -L5000  -r1 -j8 -p4
14	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_14.out -q10 -m1 -t1.5 -b100000 -L5000  -r1  -j9
15	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_15.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -jx
16	ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_16.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3
17	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_17.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -c0 
18	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_18.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 
19	ima3_stdtest -i ima3test4pop5loci.u  -o ima3_test_19.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -r5
20	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_20.out  -q10 -m1 -t1.5 -b100000 -L5000  -j01
21	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -hmg -ha0.95 -hb0.9 -jhd
21	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -hms -ha0.97 
21	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -hme -c4
22	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
23	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
24	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
25	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
26	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8 -hmg -ha0.95 -hb0.9 -jhd
27	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -hmg -ha0.95 -hb0.9 -fima3_test_28.out -jhd
28	ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_30.out  -b100000 -L5000  -j2 -c3 -g ima3_test_27.out.imapriors.txt
jh2_1	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_31.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
30	ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_32.out -q10 -m1 -t1.5 -b100000 -L5000  -j013
29	ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_33.out -q10 -m1 -t1.5 -b100000 -L5000  -j03
jh2_3	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_34.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
jh2_4	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_35.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_34.out.mcf  -jh2 
jh2_16	ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_36.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2
mpi_1	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_37.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -hmg -ha0.97 -hb0.85 -jhd
mpi_3	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_38.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn24 -hmg -ha0.97 -hb0.85 -jhd
mpi_4	mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_39.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_38.out.mcf   -hn24 -hms -ha0.97 -jhd
mpi_26	mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -hms -ha0.97 -c4
mpi_27 mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_41.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn20 -hms -ha0.97 -c4 -fima3_test_40.out 
	
17 ok
18 

ok  14,15 16, 17


4/27/2018  J. Hey

fixed some bugs revised some menu things and

	all -hf  are changed to -hm  
	
	the default heated update no longer heats the priorratio - seems to give a higher update rate and allow more heating 
	
new versions of old commands  and renumbering of commands 

ima3_stdtest -help  > ima3_test_0.out 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -hmg -ha0.95 -hb0.8 -d10 -r1 -c3  -gima3_priorfile_3pops.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_7.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hms -ha0.98 -d5
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_8.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hme -d5
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
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -hmg -ha0.95 -hb0.9
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -hms -ha0.97 
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn12 -hme -c4
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8 -hmg -ha0.95 -hb0.9
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -hmg -ha0.95 -hb0.9 -fima3_test_28.out
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_30.out  -b100000 -L5000  -j2 -c3 -g ima3_test_25.out.imapriors.txt
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_31.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_32.out -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_33.out -q10 -m1 -t1.5 -b100000 -L5000  -j03
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_34.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_35.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_34.out.mcf  -jh2 
ima3_stdtest -i ima3test3pop4loci.u  -o ima3_test_36.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_37.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -hmg -ha0.97 -hb0.85
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_38.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn24 -ha0.97 -hb0.85
mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_39.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_38.out.mcf   -hn24 -hms -ha0.97
mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_40.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn20 -hms -ha0.97 -c4
mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_41.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn20 -hms -ha0.97 -fima3_test_40.out
	

2/2/2018  J. Hey

IMa3 testbed

In windows run using ima3_stdtest.bat
in linux using ima3_stdtest.sh

both require an executable
	in windows: IMa3_stdtest.exe 
	in linux :  IMa3_stdest
		executable can be compiled using: make stdtest



There are 39 jobs, including:
	30 runs of type test_#,  where # goes from 0 thru 29

	4 runs of type test_jh2_#,  where # takes values 1,3,4 and 16  
		these runs are similar to corresponding runs test_1,test_3,test_4 and test_16
	5 runs of type test_mpi_#, where # takes values 1,3,4,26,27
		these are similar to corresponding runs test_1,test_3,test_4, test_26 and test_27
		
	set up in two batch files as of 1/17/2018

All runs produce a primary output file,  as indicated with the '-o' flag. 

In addition some runs produce additional output files (see below).

Some runs require input from previous runs (see dependencies listed below).

The following input files need to be available: 
	ima3test8pop2loci.u
	ima3test4pop2loci.u	
	ima3test3pop4loci.u
	ima3test4pop5loci.u
	pop3_nested_model_test.txt
	randNums
	ima3_priorfile_3pops.txt

The batch file is:
		ima3_stdtest.bat
		
Having weird problems copying and running these bat files
	to copy to a new dir,  open (edit) in notepad, and saveas  to new folder 

compiling and running:
	for running these tests the program must be compiled with STDTEST defined and with MPI_ENABLED defined

	In command lines listed below the program the executable has been named ima3_stdtest for clarity
	
	

test_#  runs  (all single cpu runs):
------------------------------------
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_1.out -q1 -m0.1 -t1 -b100000 -L5000 -d200 -r1 -c0 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_2.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -r1 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_3.out.mcf
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_5.out -q10 -m1 -t1.5 -r0 -vima3_test_4 -c2 -wpop3_nested_model_test.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_6.out  -b100000 -L5000 -hn10 -hfg -ha0.95 -hb0.8 -d10 -r1 -c3  -gima3_priorfile_3pops.txt
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_7.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j9 -c4 -hn50 -hfs -ha0.98 -d5
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_8.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j5 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_9.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j6 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_10.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j1 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_11.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j7 -p4
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_12.out -q10 -t1.5 -b100000 -L5000  -r1 -j8 -p4
ima3_stdtest -help  > ima3_test_13.out 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_14.out -q10 -m1 -t1.5 -b100000 -L5000  -r1  -j9
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_15.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -jx
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_16.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_17.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -c0 
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_18.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 
ima3_stdtest -i ima3test4pop5loci.u  -o ima3_test_19.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0 -r5
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_20.out  -q10 -m1 -t1.5 -b100000 -L5000  -j01
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_21.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013  -hn8 -ha0.95 -hb0.9
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_22.out  -q10 -m1 -t1.5 -b100000 -L5000  -j012
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_23.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_24.out  -q10 -m1 -t1.5 -b100000 -L5000  -j0123
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_25.out  -q10 -m1 -t1.5 -b100000 -L5000  -j23
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8 -ha0.95 -hb0.9
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -ha0.95 -hb0.9 -fima3_test_26.out
ima3_stdtest -i ima3test4pop2loci.u  -o ima3_test_28.out  -b100000 -L5000  -j2 -c3 -g ima3_test_25.out.imapriors.txt

ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_29.out  -q10 -m1 -t1.5 -b100000 -L5000  -j03
ima3_stdtest -i ima3test8pop2loci.u  -o ima3_test_30.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013


test_jh2_# runs (all single cpu runs):
--------------------------------------
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_1.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -jh2
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -jh2
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_jh2_3.out.mcf  -jh2 
ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_jh2_16.out -q10 -m1 -t1.5 -b100000 -L5000  -r1 -j3 -jh2


type test_mpi_# runs (all with 4 cpus):
---------------------------------------
C:\Program Files\Microsoft HPC Pack 2012\Bin\mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_mpi_1.out -q1 -m0.1 -t1  -b100000 -L5000 -d200 -r1 -c0 -p4 -hn12 -ha0.97 -hb0.85
C:\Program Files\Microsoft HPC Pack 2012\Bin\mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_mpi_3.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p23567 -r1245  -hn12 -ha0.97 -hb0.85
C:\Program Files\Microsoft HPC Pack 2012\Bin\mpiexec -n 4 ima3_stdtest -i ima3test3pop4loci.u -o ima3_test_mpi_4.out -q10 -m1 -t1.5 -b100000 -L5000 -d200 -c1  -p3567 -r245 -r3 -fima3_test_mpi_3.out.mcf   -hn12 -ha0.97 -hb0.85
C:\Program Files\Microsoft HPC Pack 2012\Bin\mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_mpi_26.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r7 -hn8 -ha0.95 -hb0.9
C:\Program Files\Microsoft HPC Pack 2012\Bin\mpiexec -n 4 ima3_stdtest -i ima3test4pop2loci.u -o ima3_test_mpi_27.out  -q10 -m1 -t1.5 -b100000 -L5000  -j013 -r37 -hn8 -ha0.95 -hb0.9 -fima3_test_mpi_26.out




data files:
-----------
ima3test3pop4loci.u
ima3test4pop5loci.u
ima3test8pop2loci.u
ima3test4pop2loci.u


additional required files:
--------------------------
test_5 requires pop3_nested_model_test.txt
test_6  requires  ima3_priorfile_3pops.txt


dependencies between runs:
--------------------------
test_4 must be run after test_3 has completed  (uses mcf file)
test_5 must be run after test_4 has completed (uses ti file)
test_27 must be run after test_26 has completed (uses mcf file) 
test_28 must be run after test_25 has completed  (uses  ima3_test_25.out.imapriors.txt)
test_jh2_4 must be run after test_jh2_3  has completed  (uses mcf file)
test_mpi_4 must be run after test_mpi_3 has completed   (uses mcf file)
test_mpi_27 must be run after test_mpi_26 has completed  (uses mcf file)


in addition to basic *.out files,  these other files should be produced:
------------------------------------------------------------------------

mcf files:
ima3_test_3.out.mcf.0 
ima3_test_4.out.mcf.0 
ima3_test_26.out.mcf.0
ima3_test_27.out.mcf.0

ima3_test_mpi_3.out.mcf.0  ima3_test_mpi_3.out.mcf.1  ima3_test_mpi_3.out.mcf.2  ima3_test_mpi_3.out.mcf.3  
ima3_test_mpi_4.out.mcf.0  ima3_test_mpi_4.out.mcf.1  ima3_test_mpi_4.out.mcf.2  ima3_test_mpi_4.out.mcf.3  
ima3_test_mpi_26.out.mcf.0  ima3_test_mpi_26.out.mcf.1  ima3_test_mpi_26.out.mcf.2   ima3_test_mpi_26.out.mcf.3
ima3_test_mpi_27.out.mcf.0  ima3_test_mpi_27.out.mcf.1  ima3_test_mpi_27.out.mcf.2  ima3_test_mpi_27.out.mcf.3

ima3_test_jh2_3.out.mcf.0  
ima3_test_jh2_4.out.mcf.0  


burntrend files: 
ima3_test_3.out.burntrend.out  
ima3_test_4.out.burntrend.out  

ima3_test_mpi_3.out.burntrend.out 
ima3_test_mpi_4.out.burntrend.out 

ima3_test_jh2_3.out.burntrend.out 
ima3_test_jh2_4.out.burntrend.out 

ti files:
ima3_test_4.out.ti
ima3_test_mpi_4.out.ti
ima3_test_jh2_4.out.ti


mpt files:
ima3_test_3.out.mpt
ima3_test_4.out.mpt  

ima3_test_mpi_3.out.mpt
ima3_test_mpi_4.out.mpt  

ima3_test_jh2_3.out.mpt
ima3_test_jh2_4.out.mpt  

imapriors files:
ima3_test_16.out.imapriors.txt
ima3_test_25.out.imapriors.txt
ima3_test_jh2_16.out.imapriors.txt


