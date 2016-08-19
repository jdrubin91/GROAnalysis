import node
import time
import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import load
class venn_struct:
	def __init__(self, D):
		self.D 	= D
	def get_comparison(self, *args):
		if len(args)==1:
			args 	= args[0]
			key 	= tuple((args,))
		else:
			key 	= tuple(args)
		if key in self.D:
			return self.D[key]
		return "interval_package: comparison not found in venn structure"


class interval:
	def __init__(self, l, LID=0, PID=0):
		assert len(l) > 1, "must be a list of tuples start and stop "
		self.start,self.stop,self.overlaps,self.UNIQUE,self.LID,self.PID = l[0], l[1], {}, {}, LID,PID
		if len(l)> 2:
			self.INFO 	= l[2:]
		else:
			self.INFO 	= ""
	def insert(self, i):
		if i.LID != self.LID:
			if i not in self.overlaps:
				self.overlaps[i] 	= 1
			if i.LID not in self.UNIQUE:
				self.UNIQUE[i.LID]  = list()
			if i.PID not in self.UNIQUE[i.LID]:
				self.UNIQUE[i.LID].append(i.PID)
	def get_unique(self, *args, **kwargs):
		IGNORE 	= list()
		FINDS 		= list()
		if "ignore" in kwargs:
			IGNORE 	= kwargs["ignore"]
		if len(args)==1:
			args 	= args[0]
		args 		= [a for a in args if a != self.LID ]
		
		if len(args) != len([1 for U in self.UNIQUE.keys() if U not in IGNORE]):
			return list()
		for a in args:
			if a not in self.UNIQUE:
				return list()

			for PID in self.UNIQUE[a]:
				FINDS.append((self.LID, PID))
		return FINDS

	def __str__(self):
		return str(self.start) +"-" + str(self.stop) + ":" + str(self.LID) + "_" + str(self.PID) + ":" + str(len(self.overlaps.keys()))
class overlap:
	def __init__(self,start, stop, overlaps={} ):
		self.start, self.stop, self.overlaps, self.type = start,stop, overlaps, "1"
	def insert_info(self, *args):
		for a in args:
			if a not in self.overlaps:
				self.overlaps[a] 	= 1
	def update(self, B):
		self.start 		= max(self.start, B.start)
		self.stop 		= min(self.stop, B.stop)
		if B not in self.overlaps:
			self.overlaps[B] 	= 1
	def __str__(self):
		return str(self.start) + "-" + str(self.stop) + ":" + ",".join([str(i.LID) + "_" + str(i.PID)  for i in self.overlaps.keys() ])

class comparison:
	def __init__(self, *args, **kwargs):
		self.distjoint 	= True
		self.intervals 	= list()
		self.trees 		= list() #only needed if not self.disjoint
		self.verbose 	= True
		self.error 		= False
		self.VENN_CALC 	= False
		assert len(args) >1 or len(args[0]) > 1, "needs to be at least two list of ranges to compare"
		if len(args)==1:
			args 		= args[0]
		self._build(args, kwargs)
	def _print_disjoint(self):
		if self.verbose and self.distjoint:
			print "All recieved input lists are completey disjoint"
		elif self.verbose:
			print "One or more recieved input lists overlap, constructing interval tree"
			
	def _build(self, LSTS, kwargs):
		for v in kwargs:	
			setattr(self, v, kwargs[v])
		for i,L in enumerate(LSTS):
			L.sort()
			self.intervals.append([ interval(l, LID=i, PID=j) for j,l in enumerate(L) ])
		self._check_disjoint()
		if not self.distjoint:
			for L in self.intervals:
				self.trees.append(node.tree( [ (l.start, l.stop,[l] ) for l in L ]  ))
	def _check_disjoint(self):
		l 	= 0
		while l < len(self.intervals) and self.distjoint:
			L,N,i 	= self.intervals[l], len(self.intervals[l]),1
			while i < N:
				if L[i].start < L[i-1].stop:
					self.distjoint 	= False
				if L[i].start >=L[i].stop:
					self.error 		= True
					raise TypeError, "interval_package: interval start must be less than interval stop"
				i+=1
			l+=1
		self._print_disjoint()
	def _get_overlaps_pairwaise_tree(self, A, node_B):
		overlaps 			= list()
		for a in A:
			FINDS 			= node_B.searchInterval((a.start, a.stop))
			if FINDS:
				overlaps 	   += [(max(a.start, st), min(a.stop, sp), [f for f in F]+ [a]) for st, sp, F in FINDS]
				overlaps.sort()
			
		return node.tree(overlaps)
	def _get_overlaps_pairwase(self, A, B):
		j,N 		= 0,len(B)
		overlaps 	= list()
		for a in A:
			while j < N and B[j].stop < a.start:
				j+=1
			if j < N and B[j].start < a.stop:
				O 	= None
				if hasattr(B[j], "type"):
					B[j].update(a)
					O 	= B[j]
				else:
					o_st, o_sp 	= max(a.start, B[j].start), min(a.stop, B[j].stop)
					O 			= overlap(o_st, o_sp, overlaps=dict([(a,1), (B[j],1)]))
				overlaps.append(O)				
		return overlaps
			


	def find_overlaps(self, *args, **kwargs):
		isolate 		= False
		ignore 			= list()
		if "isolate" in kwargs:
			isolate 	= kwargs["isolate"]
		if "ignore" in kwargs:
			ignore 		= kwargs["ignore"]
		if not isolate:
			assert len(args) > 1 or len(args[0]) > 1, "interval_package: need at least two lists to compare"
			if len(args)==1:
				args 	= args[0]
			searchable 	= [ self.intervals[a] for a in args if a < len(self.intervals) and a not in ignore ]
			assert len(searchable) == len(args), "interval_package: one or more of your interval comparisons is not in comparison struct"
			
			if self.distjoint:
				for i in range(len(searchable)-1):
					if i == 0:
						overlaps 	= self._get_overlaps_pairwase(searchable[i],searchable[i+1] )
					else:
						overlaps 	= self._get_overlaps_pairwase(searchable[i+1],overlaps )
				O 	= overlaps
			else:
				trees 	= [ self.trees[a] for a in args if a < len(self.intervals) and a not in ignore ]
				for i in range(len(searchable)-1):
					if i == 0:
						overlaps 	= self._get_overlaps_pairwaise_tree(searchable[i],trees[i+1] )
					else:
						overlaps 	= self._get_overlaps_pairwaise_tree(searchable[i+1], overlaps)
				O 	= overlaps.get_all()
				O.sort()
				O 	= [overlap(o[0], o[1], overlaps=dict([(I,1) for I in o[2]]))  for o in O]
			return O
		else:
			a 		= args[0]
			if not self.trees:
				self.trees 	= [node.tree( [ (l.start, l.stop,[l] ) for l in L  ]  ) for i,L in enumerate(self.intervals) if i not in ignore]
			trees 	= [ self.trees[i] for i,A in enumerate(self.intervals) if i != a  ]
			DISTINCT= list()
			for A in self.intervals[a]:
				FOUND 			= False

				for t in trees:
					if t.searchInterval((A.start, A.stop)):
						FOUND 	= True
						break
				if not FOUND:
					DISTINCT.append(A)
			return DISTINCT


	def _list_powerset(self,lst):
		# the power set of the empty set has one element, the empty set
		result = [[]]
		for x in lst:
			result.extend([subset + [x] for subset in result])
		return result
	def get_counts(self,*args, **kwargs):
		IGNORE 	= list()
		if "ignore" in kwargs:
			IGNORE 	= kwargs["ignore"]
		if len(args)==1:
			args 	= args[0]
		C 	= dict([(a,0) for a in args])
		for a in args:
			assert a < len(self.intervals), "asked for the wrong list IDS must be less than " + str(self.intervals)
			for I in self.intervals[a]:
				FINDS 	= I.get_unique(args, ignore=IGNORE)
				C[a]+=int(len(FINDS) > 0 )
		return C
	def get_isolated(self, a, ignore=list()):
		assert a < len(self.intervals), "asked for the wrong list IDS must be less than " + str(self.intervals)
		if self.VENN_CALC:
			CT = list()
			for I in self.intervals[a]:
				FINDS 	= len([1 for II in I.UNIQUE.keys() if II not in ignore ] )
				if FINDS == 0:
					CT.append(I)
		else:
			return self.find_overlaps(a, isolate=True, ignore=ignore)
		return CT
	def _draw_venn_dendrogram(self, *args, **kwargs):
		args 	= args[0]
		if "labels" not in kwargs:
			labels 	= ["set " + str(a+1) for a in args]
		else:
			labels 	= kwargs["labels"]	
		if len(args)==3:
			F 	= plt.figure(figsize=(5,5))
			ax1 = F.add_subplot(1,1,1)
			ax1.set_title("Venn Diagram of Overlapping Events")
			
			circle1=plt.Circle((0.5,1.0/3.),.2,facecolor='r', alpha=0.333, linewidth=1, edgecolor="black",linestyle="dashed" )
			circle1.set_radius(0.3)

			circle2=plt.Circle((1./3.,2/3.),.2,facecolor='b', alpha=0.333, linewidth=1, edgecolor="black", linestyle="dashed" )
			circle2.set_radius(0.3)
			circle3=plt.Circle((2./3.,2/3.),.2,facecolor='g', alpha=0.333, linewidth=1, edgecolor="black" ,linestyle="dashed" )
			circle3.set_radius(0.3)
			i,j 	= 0,0
			IGNORE 	= [i for i in range(len(self.intervals)) if i not in args  ] 
			tri  	= self.get_counts(args[0],args[1],args[2], ignore=IGNORE) 
			pair_0_1 = self.get_counts(args[0],args[1], ignore=IGNORE )
			pair_0_2 = self.get_counts(args[0],args[2], ignore=IGNORE )
			pair_1_2 = self.get_counts(args[1],args[2], ignore=IGNORE )
			ax1.text(0.49,0.48,  "\n".join([str(x) for x in tri.values() ]) , fontsize=10, color="w")
			ax1.text(0.64,0.43,  "\n".join([str(x) for x in pair_0_2.values() ]), fontsize=12, color="w")
			ax1.text(0.3,0.43,  "\n".join([str(x) for x in pair_0_1.values() ]), fontsize=12, color="w")
			ax1.text(0.49,2/3.,  "\n".join([str(x) for x in pair_1_2.values() ]), fontsize=12, color="w")
			ax1.text(0.49,0.25,  str( len(self.get_isolated(args[0], ignore=IGNORE) )), fontsize=18, color="w")
			ax1.text(0.25,2/3.,  str( len(self.get_isolated(args[1], ignore=IGNORE) )), fontsize=18, color="w")
			ax1.text(0.75,2/3.,  str( len(self.get_isolated(args[2], ignore=IGNORE)  )), fontsize=18, color="w")
			ax1.legend([circle1, circle2, circle3], labels, loc=(0.8,0.1), fontsize=12)
			F.gca().add_artist(circle1)
			F.gca().add_artist(circle2)
			F.gca().add_artist(circle3)

		else:
			IGNORE 	= [i for i in range(len(self.intervals)) if i not in args  ] 
			F 	= plt.figure(figsize=(5,5))
			ax1 = F.add_subplot(1,1,1)
			ax1.set_title("Venn Diagram of Overlapping Events")
			
			circle1=plt.Circle((1/3.,0.5),.2,facecolor='r', alpha=0.333, linewidth=1, edgecolor="black",linestyle="dashed" )
			circle1.set_radius(0.3)

			circle2=plt.Circle((2/3., 0.5),.2,facecolor='b', alpha=0.333, linewidth=1, edgecolor="black",linestyle="dashed" )
			circle2.set_radius(0.3)
			pair_0_1 = self.get_counts(args[0],args[1], ignore=IGNORE )
			ax1.text(0.5,0.475,  "\n".join([str(x) for x in pair_0_1.values() ]), fontsize=12, color="w")
			ax1.text(0.25,0.5,  str( len(self.get_isolated(args[0], ignore=IGNORE) )), fontsize=18, color="w")
			ax1.text(0.75,0.5,  str( len(self.get_isolated(args[1], ignore=IGNORE) )), fontsize=18, color="w")
			F.gca().add_artist(circle1)
			F.gca().add_artist(circle2)
			ax1.legend([circle1, circle2], labels, loc=(0.8,0.1), fontsize=12)
			
		ax1.set_xticklabels([])
		ax1.set_yticklabels([])


		plt.show()
	

	def compute_venn(self, *args, **kwargs):
		if len(args)==1:
			args 	= args[0]
		display 	= False
		if "display" in kwargs:
			display 	= kwargs['display']
		power_set 		= [ i for i in self._list_powerset( args ) if (len(i) > 1)]
		A 				= {}
		self.VENN_CALC 	= True
		for c in power_set:
			
			O 		= self.find_overlaps(c)
			for o in O:
				for i in  o.overlaps :
					for j in o.overlaps:
						if i != j:
							i.insert(j)
							j.insert(i)
			A[tuple(c)] = O
		power_set 	= [(len(i),tuple(i) ) for i in power_set]
		power_set.sort()
		power_set[::-1]
		power_set 	= [c for i,c in power_set]
		for a in args:
			A[tuple((a,))] = self.get_isolated(a)

		if display:
			assert 1< len(args) < 4,"interval_package: display only supported for 2 or 3 comparisons"
			if "labels" in kwargs:
				self._draw_venn_dendrogram(args, labels=kwargs["labels"])
			else:
				self._draw_venn_dendrogram(args)
				
		return venn_struct(A)


	



