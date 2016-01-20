import intervals,load
#===================================================================
#TOY TEST DATA
#needs to be a list of tuples where the first two elements are the 
#start and stop of the intervals. The rest of the tuple can be there
#or not and can be any size

A 	= [(1,5, "A1","hi"), (4,10, "A2"),  (13,15, "A3"), (32, 34, "A4"), (61,68, "A5")]
B 	= [(1,6), (7,15, "B2"),  (16,17, "B3" ), (62,69, "B4") ]
C 	= [(2,6, "C1"), (18,20, "C2"),  (21,23, "C3"), (25, 29), (31, 35)]
D 	= [(2,7, "D1"), (12, 17, "D2"), (61,65, "D3")]

#first thing is the initialize the intervals data structure
#this will figure out if an interval tree is the best
#algorithm or just a regular list comparison
ST 	= intervals.comparison((A,B,C, D) )

#===================================================================
#there are three important methods that you can use
#FIRST ask for overlaps between any combination of the lists

OVERLAPS_0_1 	= ST.find_overlaps(0,1)
OVERLAPS_0_1_2 	= ST.find_overlaps(0,1,2)

#===================================================================
#now ST.find_overlaps will return a list of overlap classes
#the overlap class will have a start and stop location of the 
#overlap event (acessed by overlap.start and overlap.stop)
#the overlap class will also have a attribute called overlaps
#with is a dictionary where the key is the original intervals
#that overlap it. These intervals have the original start and stop
#attributes and the info e.g. ->
print "---------------------------------------------"
print "There were", len(OVERLAPS_0_1_2), "overlaps between lists 0,1,2" 
for O in OVERLAPS_0_1_2:
	print "Overlap Instance: ", O, O.start, O.stop, len(O.overlaps.keys())
	print "Intervals within: "
	for interval_original in O.overlaps:
	        print interval_original.INFO == ''
		print interval_original, interval_original.start, interval_original.stop, interval_original.INFO
print "---------------------------------------------"

#===================================================================
#SECOND you can also get those intervals that don't overlap anybody
#between the lists used to construct ST, use ST.get_isolated
UNIQUE_2 	= ST.get_isolated(2)

print "---------------------------------------------"
print "There were", len(UNIQUE_2), "unique elements in LST 3" 
for I in UNIQUE_2:
	print I, I.start, I.stop, I.INFO
print "---------------------------------------------"

#===================================================================
#THIRD you can compute all possible, mutually exclusive overlaps
#i.e. a Venn Diagram

V 									= ST.compute_venn(0,1, 2,3, display=False )
UNIQUE_2 							= V.get_comparison(2)
ALL_NOT_MUTUAL_EXCLUSIVE_OVERLAPS 	= V.get_comparison(0,1,3)

#IMPORTANT, althrough the calculations in the displayed venn diagram are correct 
#(they can be asked by the get_counts() method see below)
#the data structure returns all overlaps and doesnt remove duplicates

MUTUAL_EXCLUSE_COUNTS 				= ST.get_counts(0,1) #not ignoring lists 2 and 3
print "---------------------------------------------"
print "Mutual exclusive counts not ignoring 2 and 3"
print MUTUAL_EXCLUSE_COUNTS

MUTUAL_EXCLUSE_COUNTS 				= ST.get_counts(0,1,ignore=(2,3)) #not ignoring lists 2 and 3
print "---------------------------------------------"
print "Mutual exclusive counts ignoring 2 and 3"
print MUTUAL_EXCLUSE_COUNTS

print "---------------------------------------------"
print "There were", len(ALL_NOT_MUTUAL_EXCLUSIVE_OVERLAPS), "mutually exclusive overlaps between lists 0,1,3" 
for O in ALL_NOT_MUTUAL_EXCLUSIVE_OVERLAPS:
	print "Overlap Instance: ", O, O.start, O.stop, len(O.overlaps.keys())
	print "Intervals within: "
	for interval_original in O.overlaps:
		print interval_original, interval_original.start, interval_original.stop, interval_original.INFO
print "---------------------------------------------"

#===================================================================
#I also wrote a plotter to show the venn diagram that depends 
#only matplotlib
V 							= ST.compute_venn(0,1, 2,  display=True , labels=("A", "B", "C"))
V 							= ST.compute_venn(0,  2,  display=True, labels=("A", "B", "C") )
#===================================================================
#go ahead and test on your on bed files

REAL_INTERVALS 				= False
if REAL_INTERVALS:
	BED_FILE_DIRECTORY 	= "/Users/joazofeifa/Lab/gro_seq_files/"
	BED_FILE_1 			= "Allen2014/EMG_out_files/Allen2014_DMSO2_3-4_bidirectional_hits_intervals.bed"
	BED_FILE_2 			= "Andersson2014/EMG_out_files/SRR1596501_Andersson2014-4_bidirectional_hits_intervals.bed"
	BED_FILE_3 			= "Luo2014/EMG_out_files/SRR1015583_Luo2014-3_bidirectional_hits_intervals.bed"
	REFSEQ_FILE 		= "/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt"


	R 				= load.refseq_table("/Users/joazofeifa/Lab/genome_files/RefSeqHG19.txt")
	#R 				= None
	Allen2014 		= load.bed_file("/Users/joazofeifa/Lab/gro_seq_files/Allen2014/EMG_out_files/Allen2014_DMSO2_3-4_bidirectional_hits_intervals.bed",FILTER=R)["chr1"]
	Anderson2014 	= load.bed_file("/Users/joazofeifa/Lab/gro_seq_files/Andersson2014/EMG_out_files/SRR1596501_Andersson2014-4_bidirectional_hits_intervals.bed",FILTER=R)["chr1"]
	Luo2014 		= load.bed_file("/Users/joazofeifa/Lab/gro_seq_files/Luo2014/EMG_out_files/SRR1015583_Luo2014-3_bidirectional_hits_intervals.bed",FILTER=R)["chr1"]


	ST 	= intervals.comparison((Allen2014,Anderson2014,Luo2014), verbose=False)
	ST.compute_venn(0,1,2, display=True )
	ST.compute_venn(0,1 , display=True )

