#####
Tests
#####

Below we provide example output for many of the tests that can be run in the SmallK library. This will provide guidance for identifying issues with your installation. Of course results will vary on different machines.

*******************
SmallK Test Results
*******************

After building the smallk library, the `make check` command will run a bash script that performs a series of tests on the code.  This is a sample output of those tests:

.. code-block:: cpp

	Build configuration: release
	sh tests/scripts/test_smallk.sh ../xdata_github/smallk_data/ | tee smallk_test_results.txt 
	*****************************************************
	*                                                   *
	*           Testing the smallk interface.           *
	*                                                   *
	*****************************************************
	WARNING: Could not achieve THREAD_MULTIPLE support.
	Smallk major version: 1
	Smallk minor version: 6
	Smallk patch level:   0
	Smallk version string: 1.6.0
	Loading matrix...

	Running NMF-BPP...

	Initializing matrix W...
	Initializing matrix H...

	                parameters: 

		         algorithm: Nonnegative Least Squares with Block Principal Pivoting
		stopping criterion: Ratio of Projected Gradients
		            height: 12411
		             width: 7984
		                 k: 8
		           miniter: 1
		           maxiter: 5000
		               tol: 0.005
		        matrixfile: ../xdata_github/smallk_data/reuters.mtx
		        maxthreads: 8

	1:	progress metric: 	(min_iter)
	2:	progress metric:	0.35826
	3:	progress metric:	0.172127
	4:	progress metric:	0.106297
	5:	progress metric:	0.0696424
	6:	progress metric:	0.0538889
	7:	progress metric:	0.0559478
	8:	progress metric:	0.0686117
	9:	progress metric:	0.0788641
	10:	progress metric:	0.0711522
	20:	progress metric:	0.00568349

	Solution converged after 22 iterations.

	Elapsed wall clock time: 0.633 sec.

	Writing output files...

	Running HierNmf2...

	Loading dictionary...

		        parameters: 

		            height: 12411
		             width: 7984
		        matrixfile: ../xdata_github/smallk_data/reuters.mtx
		          dictfile: ../xdata_github/smallk_data/reuters_dictionary.txt
		               tol: 0.0001
		           miniter: 1
		           maxiter: 5000
		          maxterms: 5
		        maxthreads: 8
	[1] [2] [3] [4] 

	Elapsed wall clock time: 551 ms.
	9/9 factorizations converged.

	Writing output files...
	W matrix test passed
	H matrix test passed
	*****************************************************
	*                                                   *
	*            Testing the preprocessor.              *
	*                                                   *
	*****************************************************

	      Command line options: 

		             indir: ../xdata_github/smallk_data/
		            outdir: current directory
		     docs_per_term: 3
		     terms_per_doc: 5
		          max_iter: 1000
		         precision: 4
		      boolean_mode: 0

	Loading input matrix ../xdata_github/smallk_data/matrix.mtx
		Input file load time: 1.421s.

	Starting iterations...
		[1] height: 39771, width: 11237, nonzeros: 877453
	Iterations finished.
		New height: 39727
		New width: 11237
		New nonzero count: 877374
	Processing time: 0.063s.

	Writing output matrix 'reduced_matrix.mtx'
	Output file write time: 2.189s.
	Writing dictionary file 'reduced_dictionary.txt'
	Writing documents file 'reduced_documents.txt'
	Dictionary + documents write time: 0.083s.
	preprocessor matrix test passed
	preprocessor dictionary test passed
	preprocessor documents test passed
	*****************************************************
	*                                                   *
	*            Testing the NMF routines.              *
	*                                                   *
	*****************************************************
	WARNING: Could not achieve THREAD_MULTIPLE support.
	Loading matrix...
	Initializing matrix W...
	Initializing matrix H...

	      Command line options: 

		         algorithm: Nonnegative Least Squares with Block Principal Pivoting
		stopping criterion: Ratio of Projected Gradients
		            height: 12411
		             width: 7984
		                 k: 8
		           miniter: 1
		           maxiter: 5000
		               tol: 0.005
		          tolcount: 1
		           verbose: 1
		         normalize: 1
		      outprecision: 6
		        matrixfile: ../xdata_github/smallk_data//reuters.mtx
		          infile_W: ../xdata_github/smallk_data//nmf_init_w.csv
		          infile_H: ../xdata_github/smallk_data//nmf_init_h.csv
		         outfile_W: w.csv
		         outfile_H: h.csv
		        maxthreads: 8

	1:	progress metric: 	(min_iter)
	2:	progress metric:	0.35826
	3:	progress metric:	0.172127
	4:	progress metric:	0.106297
	5:	progress metric:	0.0696424
	6:	progress metric:	0.0538889
	7:	progress metric:	0.0559478
	8:	progress metric:	0.0686117
	9:	progress metric:	0.0788641
	10:	progress metric:	0.0711522
	20:	progress metric:	0.00568349

	Solution converged after 22 iterations.

	Elapsed wall clock time: 0.673 sec.

	Writing output files...
	NMF W matrix test passed
	NMF H matrix test passed
	*****************************************************
	*                                                   *
	*                Testing hierclust.                 *
	*                                                   *
	*****************************************************
	-------------------------------
	                               
	  Reuters matrix, 12 clusters  
	                               
	-------------------------------
	WARNING: Could not achieve THREAD_MULTIPLE support.
	loading dictionary...
	loading matrix...

	     Command line options: 

		            height: 12411
		             width: 7984
		        matrixfile: ../xdata_github/smallk_data//reuters.mtx
		           initdir: ../xdata_github/smallk_data//test/matrices.reuters/
		          dictfile: ../xdata_github/smallk_data//reuters_dictionary.txt
		        assignfile: assignments_12.csv
		            format: XML
		          treefile: tree_12.xml
		          clusters: 12
		               tol: 0.0001
		            outdir: 
		           miniter: 1
		           maxiter: 5000
		          maxterms: 5
		        maxthreads: 8
		        unbalanced: 0.1
		   trial_allowance: 3
		              flat: 0
		           verbose: 1

	[1] [2] [3] [4] [5] [6] dropping 20 items ...
	[7] [8] [9] [10] [11] 

	Elapsed wall clock time: 2.758 s.
	26/26 factorizations converged.

	Writing output files...
	XML file test passed
	assignment file test passed
	-------------------------------
	                               
	  20News matrix, 15 clusters  
	                               
	-------------------------------
	WARNING: Could not achieve THREAD_MULTIPLE support.
	loading dictionary...
	loading matrix...

	     Command line options: 

		            height: 39727
		             width: 11237
		        matrixfile: ../xdata_github/smallk_data//news20.mtx
		           initdir: ../xdata_github/smallk_data//test/matrices.20news/
		          dictfile: ../xdata_github/smallk_data//news20_dictionary.txt
		        assignfile: assignments_15.csv
		            format: XML
		          treefile: tree_15.xml
		          clusters: 15
		               tol: 0.0001
		            outdir: 
		           miniter: 1
		           maxiter: 5000
		          maxterms: 5
		        maxthreads: 8
		        unbalanced: 0.1
		   trial_allowance: 3
		              flat: 0
		           verbose: 1

	[1] [2] [3] dropping 30 items ...
	[4] [5] dropping 132 items ...
	[6] [7] [8] [9] [10] dropping 41 items ...
	[11] dropping 51 items ...
	[12] dropping 22 items ...
	[13] dropping 85 items ...
	[14] 

	Elapsed wall clock time: 10.308 s.
	41/41 factorizations converged.

	Writing output files...
	XML file test passed
	assignment file test passed
	*****************************************************
	*                                                   *
	*              Testing flatclust.                   *
	*                                                   *
	*****************************************************
	WARNING: Could not achieve THREAD_MULTIPLE support.
	loading dictionary...
	loading matrix...
	Initializing matrix W...
	Initializing matrix H...

	     Command line options: 

		            height: 256
		             width: 256
		        matrixfile: ../xdata_github/smallk_data//rnd_256_256.csv
		          infile_W: ../xdata_github/smallk_data//flatclust_init_w.csv
		          infile_H: ../xdata_github/smallk_data//flatclust_init_h.csv
		          dictfile: ../xdata_github/smallk_data//reuters_dictionary.txt
		        assignfile: assignments_16.csv
		         fuzzyfile: assignments_fuzzy_16.csv
		            format: XML
		         clustfile: clusters_16.xml
		         algorithm: HALS
		          clusters: 16
		               tol: 0.0001
		            outdir: 
		           miniter: 1
		           maxiter: 5000
		          maxterms: 5
		        maxthreads: 8
		           verbose: 1

	1:	progress metric: 	(min_iter)
	2:	progress metric:	0.635556
	3:	progress metric:	0.490817
	4:	progress metric:	0.479135
	5:	progress metric:	0.474986
	6:	progress metric:	0.44968
	7:	progress metric:	0.422542
	8:	progress metric:	0.407662
	9:	progress metric:	0.395145
	10:	progress metric:	0.379238
	20:	progress metric:	0.272868
	30:	progress metric:	0.168386
	40:	progress metric:	0.109147
	50:	progress metric:	0.0767327
	60:	progress metric:	0.0488545
	70:	progress metric:	0.036226
	80:	progress metric:	0.0307648
	90:	progress metric:	0.0266116
	100:	progress metric:	0.0226963
	110:	progress metric:	0.0188616
	120:	progress metric:	0.0158307
	130:	progress metric:	0.0137605
	140:	progress metric:	0.0127888
	150:	progress metric:	0.0123962
	160:	progress metric:	0.0124734
	170:	progress metric:	0.0123563
	180:	progress metric:	0.0122163
	190:	progress metric:	0.0120643
	200:	progress metric:	0.0117647
	210:	progress metric:	0.0114894
	220:	progress metric:	0.0110467
	230:	progress metric:	0.0107816
	240:	progress metric:	0.0105239
	250:	progress metric:	0.0103824
	260:	progress metric:	0.0100915
	270:	progress metric:	0.00965073
	280:	progress metric:	0.00938526
	290:	progress metric:	0.00914129
	300:	progress metric:	0.00896701
	310:	progress metric:	0.00886729
	320:	progress metric:	0.00841059
	330:	progress metric:	0.007793
	340:	progress metric:	0.00740095
	350:	progress metric:	0.00708869
	360:	progress metric:	0.00683069
	370:	progress metric:	0.00672093
	380:	progress metric:	0.00687906
	390:	progress metric:	0.00703777
	400:	progress metric:	0.00721928
	410:	progress metric:	0.00729384
	420:	progress metric:	0.00718332
	430:	progress metric:	0.00722893
	440:	progress metric:	0.00726766
	450:	progress metric:	0.00739665
	460:	progress metric:	0.00769819
	470:	progress metric:	0.00814673
	480:	progress metric:	0.008566
	490:	progress metric:	0.00877955
	500:	progress metric:	0.00884221
	510:	progress metric:	0.0088057
	520:	progress metric:	0.00852345
	530:	progress metric:	0.00797952
	540:	progress metric:	0.00749354
	550:	progress metric:	0.00689316
	560:	progress metric:	0.00623287
	570:	progress metric:	0.00576619
	580:	progress metric:	0.00541125
	590:	progress metric:	0.00501715
	600:	progress metric:	0.00466547
	610:	progress metric:	0.00432811
	620:	progress metric:	0.00412669
	630:	progress metric:	0.00383406
	640:	progress metric:	0.00352802
	650:	progress metric:	0.00331556
	660:	progress metric:	0.00315735
	670:	progress metric:	0.00304253
	680:	progress metric:	0.00296627
	690:	progress metric:	0.00289013
	700:	progress metric:	0.00279647
	710:	progress metric:	0.00271036
	720:	progress metric:	0.00261087
	730:	progress metric:	0.0025158
	740:	progress metric:	0.00245123
	750:	progress metric:	0.00237435
	760:	progress metric:	0.00231126
	770:	progress metric:	0.00228199
	780:	progress metric:	0.00227623
	790:	progress metric:	0.00228185
	800:	progress metric:	0.00227993
	810:	progress metric:	0.00228216
	820:	progress metric:	0.00228018
	830:	progress metric:	0.00229096
	840:	progress metric:	0.00232403
	850:	progress metric:	0.00234957
	860:	progress metric:	0.00227868
	870:	progress metric:	0.00210786
	880:	progress metric:	0.00195462
	890:	progress metric:	0.00183587
	900:	progress metric:	0.00173358
	910:	progress metric:	0.0016405
	920:	progress metric:	0.00156422
	930:	progress metric:	0.00150835
	940:	progress metric:	0.00146594
	950:	progress metric:	0.00143261
	960:	progress metric:	0.00137378
	970:	progress metric:	0.00131989
	980:	progress metric:	0.00126626
	990:	progress metric:	0.0012164
	1000:	progress metric:	0.00117061
	1010:	progress metric:	0.00112539
	1020:	progress metric:	0.00108626
	1030:	progress metric:	0.00105192
	1040:	progress metric:	0.00102131
	1050:	progress metric:	0.000992069
	1060:	progress metric:	0.000965259
	1070:	progress metric:	0.000938949
	1080:	progress metric:	0.000911962
	1090:	progress metric:	0.000884505
	1100:	progress metric:	0.000854904
	1110:	progress metric:	0.000820121
	1120:	progress metric:	0.000785245
	1130:	progress metric:	0.000752513
	1140:	progress metric:	0.000723279
	1150:	progress metric:	0.000697698
	1160:	progress metric:	0.000680904
	1170:	progress metric:	0.000652152
	1180:	progress metric:	0.000628268
	1190:	progress metric:	0.000612413
	1200:	progress metric:	0.000596834
	1210:	progress metric:	0.000580674
	1220:	progress metric:	0.000556549
	1230:	progress metric:	0.000535666
	1240:	progress metric:	0.00051492
	1250:	progress metric:	0.000496234
	1260:	progress metric:	0.000481147
	1270:	progress metric:	0.000461294
	1280:	progress metric:	0.000440802
	1290:	progress metric:	0.000419049
	1300:	progress metric:	0.000398007
	1310:	progress metric:	0.000376203
	1320:	progress metric:	0.000355811
	1330:	progress metric:	0.00033729
	1340:	progress metric:	0.000318932
	1350:	progress metric:	0.000302528
	1360:	progress metric:	0.000287961
	1370:	progress metric:	0.00027486
	1380:	progress metric:	0.00026403
	1390:	progress metric:	0.000255504
	1400:	progress metric:	0.000248646
	1410:	progress metric:	0.000242996
	1420:	progress metric:	0.000239243
	1430:	progress metric:	0.000236852
	1440:	progress metric:	0.000235313
	1450:	progress metric:	0.000234465
	1460:	progress metric:	0.000234154
	1470:	progress metric:	0.000234253
	1480:	progress metric:	0.00023487
	1490:	progress metric:	0.000237223
	1500:	progress metric:	0.000240043
	1510:	progress metric:	0.000243896
	1520:	progress metric:	0.00024867
	1530:	progress metric:	0.000253981
	1540:	progress metric:	0.000260239
	1550:	progress metric:	0.000266795
	1560:	progress metric:	0.000273529
	1570:	progress metric:	0.000280678
	1580:	progress metric:	0.000287273
	1590:	progress metric:	0.000292288
	1600:	progress metric:	0.000296475
	1610:	progress metric:	0.000299556
	1620:	progress metric:	0.00030244
	1630:	progress metric:	0.000306148
	1640:	progress metric:	0.000310299
	1650:	progress metric:	0.000314674
	1660:	progress metric:	0.000319052
	1670:	progress metric:	0.000323906
	1680:	progress metric:	0.000329536
	1690:	progress metric:	0.000335913
	1700:	progress metric:	0.000342834
	1710:	progress metric:	0.000351167
	1720:	progress metric:	0.000352515
	1730:	progress metric:	0.000348749
	1740:	progress metric:	0.000345684
	1750:	progress metric:	0.000343139
	1760:	progress metric:	0.000340867
	1770:	progress metric:	0.000339052
	1780:	progress metric:	0.000337038
	1790:	progress metric:	0.000335244
	1800:	progress metric:	0.000333452
	1810:	progress metric:	0.000332111
	1820:	progress metric:	0.000330198
	1830:	progress metric:	0.000325983
	1840:	progress metric:	0.000321473
	1850:	progress metric:	0.000316999
	1860:	progress metric:	0.000312054
	1870:	progress metric:	0.000305176
	1880:	progress metric:	0.000294684
	1890:	progress metric:	0.000284482
	1900:	progress metric:	0.000274905
	1910:	progress metric:	0.000265684
	1920:	progress metric:	0.000256761
	1930:	progress metric:	0.000248203
	1940:	progress metric:	0.000239613
	1950:	progress metric:	0.000230677
	1960:	progress metric:	0.00022218
	1970:	progress metric:	0.000214089
	1980:	progress metric:	0.00020621
	1990:	progress metric:	0.000196915
	2000:	progress metric:	0.000187712
	2010:	progress metric:	0.000179199
	2020:	progress metric:	0.00017137
	2030:	progress metric:	0.000164158
	2040:	progress metric:	0.000157751
	2050:	progress metric:	0.000152485
	2060:	progress metric:	0.000147217
	2070:	progress metric:	0.000142083
	2080:	progress metric:	0.000137148
	2090:	progress metric:	0.000132379
	2100:	progress metric:	0.000127922
	2110:	progress metric:	0.000123617
	2120:	progress metric:	0.000119548
	2130:	progress metric:	0.000115684
	2140:	progress metric:	0.000111997
	2150:	progress metric:	0.000108389
	2160:	progress metric:	0.000104838
	2170:	progress metric:	0.000101387

	Solution converged after 2175 iterations.

	Elapsed wall clock time: 1.022 sec.

	XML file test passed
	assignment file test passed
	fuzzy assignment file test passed
	***** SmallK: All tests passed. *****

