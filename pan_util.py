#!/usr/bin/python
import sys
from lockgraph import *

if_print = 1
throu1 = []
throu3 = []
throu_diff = []
classList = []
n_item1 = []
n_item3 = []
item_diff = []
#def update_classL(threadL_todelete,classL):
#	for c in classL:
def keep_ts_in_range(tryDict,acqDict,relDict,start_ts,end_ts):
	newRelDict = {}
	newTryDict = {}
	newAcqDict = {}
	for k,v in relDict.iteritems(): #v is a list of the tuple (lockid,timstamp)
		newRelDict[k] = []
		newTryDict[k] = []
		newAcqDict[k] = []
		i = 0
		for t in v:
			if t[1] >= start_ts and t[1] <= end_ts: #the tuple index starts from 0, if the timestamp falls in the rage
				newRelDict[k].append(t)
				newTryDict[k].append(tryDict[k][i])
				newAcqDict[k].append(acqDict[k][i])
			i = i + 1	
		#print "in tryDict, original size of list for key ", k, len(tryDict[k])
		#print "in tryDict, new size of list for key ", k, len(newTryDict[k])
		t = 0
		for x in newTryDict[k]:
			if x[0] == 0:
				t = t+1
		#print "number of entries for lock 0 in key", k, ": ",t
		if newTryDict[k] == [] or newAcqDict[k] == [] or newRelDict[k] == []:
			del newTryDict[k]
			del newAcqDict[k]
			del newRelDict[k]
	return (newTryDict,newAcqDict,newRelDict)
		
def minT(tList,index): #tList is a list of list
	minT = (tList[0])[0]
	for t in tList:
		for s in t:
			if s[index] < minT[index]:
				minT = s
	return minT

def maxT(tList,index): #tList is a list of list
	maxT = (tList[0])[0]
	for t in tList:
		for s in t:
			if s[index] > maxT[index]:
				maxT = s
	return maxT

def getClassL (newThreadL,originalThreadL,n):
	stage2 = []
	stage3 = []
	stage4 = []
	newClassL = map(originalThreadL.index,newThreadL)
	j = 0
	for i in newClassL:
		if i < n:
			stage2.append(j)
		elif i < 2*n and i >=n:
			stage3.append(j)
		elif i >= 2*n and i < 3*n:
			stage4.append(j)
		j = j+1
	return [stage2,stage3,stage4]

def count_entries(tryDict, classID, lockID):
	#print "in count_entries", classList
	count = 0
	for k,v in tryDict.iteritems(): #v is a list of the tuple (lockid,timstamp)
		if k in classList[classID]:
			for x in v:
				if x[0] == lockID:
					count = count + 1
	return count 
	
def slotted_analyzed(start,end,nDiv,slot_size,tryDic,acqDic,relDic,namesD,n,if_class1_finished):
	if nDiv:
		slot_size = (end - start)/nDiv
	else:
		nDiv = (end-start)/slot_size
	m1 = [] #class 1 measured
	m3 = [] #class 3 measured
	a1 = [] #class 1 analysis
	a3 = [] #class 3 analysis
	s1 = [] #class 1 service time
	s3 = [] #class 3 service time
	t = []

	for i in range(nDiv):
		print "slot ", i, ":"
		slot_start = start + i*slot_size 
		if i == nDiv - 1:
			slot_end = end
		else:
			slot_end = slot_start + slot_size

		#print "slot starts: ",slot_start, " ends: ",slot_end, " slot size:", slot_end - slot_start
		(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,slot_start,slot_end)
		newClassL = getClassL (newTryDict.keys(),tryDic.keys(),n)
		#print "newClassL:", newClassL

		#get the number of lock accesses for a specific class for a specific lock
		print "count class 0:",count_entries(newTryDict,0,0)
		print "count class 1:",count_entries(newTryDict,2,0)

		res = multi_analyze(newTryDict, newAcqDict, newRelDict, namesD, newClassL)
		measured = res[4][0] #first lock waiting time
		#print "measured:", measured
		newnewClassL = filter(lambda x: x != [], newClassL)
		#print "newClassL:", newClassL

		serv = multi_serviceTime(newTryDict,newAcqDict,newRelDict,namesD,newnewClassL,0)
		#print "newnewClassL:", newnewClassL
		classVector = tuple(map(len,newnewClassL))
		#analyzed = mva_multiclass(*res[:4])[0]
		#print "classVector:", classVector
		ana = mva_multiclass(res[0],res[1],classVector,res[3])
		analyzed = ana[0]
		#print "Multiclass MVA:", analyzed[1]
		#print "printing service times:"
		#print "analyzed[2][1]:",ana[2][1]

		n1 = count_entries(newTryDict,0,0)
		n_item1.append(n1)
		n3 = count_entries(newTryDict,2,0)
		n_item3.append(n3)
		item_diff.append(sum(n_item1) - sum(n_item3))

		
		if len(measured) == 3:
			throu1.append(ana[2][1][0]*1e9)
			throu3.append(ana[2][1][2]*1e9)
			throu_diff.append(sum(throu1) - sum(throu3))

			#print "throughput 1: ",ana[2][1][0]*1e9,"throughput 3: ",ana[2][1][2]*1e9,"difference: ",(ana[2][1][0] - ana[2][1][2])*1e9
		if len(measured) == 2:
			if not if_class1_finished:
				throu1.append(ana[2][1][0]*1e9)
			throu3.append(ana[2][1][1]*1e9)
			#throu_diff.append(ana[2][1][0]*1e9 - ana[2][1][1]*1e9)
			throu_diff.append(sum(throu1) - sum(throu3))
			#print "throughput 1: ",ana[2][1][0]*1e9,"throughput 3:",ana[2][1][1]*1e9,"throughput difference",(ana[2][1][0] - ana[2][1][1])*1e9
		if len(measured) == 1:
			throu3.append(ana[2][1][0]*1e9)
			#throu_diff.append(-ana[2][1][0]*1e9)
			throu_diff.append(sum(throu1)-sum(throu3))
			#print "throughput of consumer",ana[2][1][0]*1e9

		if len(measured) == 3:
			m1.append(measured[0])
			m3.append(measured[2])
			a1.append(analyzed[1][0])
			a3.append(analyzed[1][2])
			s1.append(serv[0][0])
			s3.append(serv[2][0])
		elif len(measured) == 2:
			m1.append(measured[0])
			m3.append(measured[1])
			a1.append(analyzed[1][0])
			a3.append(analyzed[1][1])
			s1.append(serv[0][0])
			s3.append(serv[1][0])
		elif len(measured) == 1:
			m3.append(measured[0])
			a3.append(analyzed[1][0])
			s3.append(serv[0][0])
	if if_print:
		print "m1:"
		for i in m1:
			print i 
		print "m3:"
		for i in m3:
			print i 
		print "a1:"
		for i in a1:
			print i 
		print "a3:"
		for i in a3:
			print i
		print "s1:"
		for i in s1:
			print i 
		print "s3:"
		for i in s3:
			print i
	return slot_size

def get_derivative_analysis(n):
	#(tryDic, acqDic, relDic, namesD) = parseDicTuple ('./dedup_run_on_halvan/dedup_10th_32c_110805/dedup.native.10th.halvan')
	print "Parsing input file:",sys.argv[1]
	#(tryDic, acqDic, relDic, namesD) = parseDicTuple ('./dedup_run_on_halvan/dedup_10th_32c_random_input/dedup.native.10th.halvan.random')
	(tryDic, acqDic, relDic, namesD) = parseDicTuple (sys.argv[1])
	map(lambda x:x[0],tryDic.values())
	#print "tryDic keys:", tryDic.keys()
	#classList = [tryDic.keys()[0:n],tryDic.keys()[n:2*n],tryDic.keys()[2*n:3*n]] #classList has the thread ids
	classList.append(tryDic.keys()[0:n])
	classList.append(tryDic.keys()[n:2*n])
	classList.append(tryDic.keys()[2*n:3*n]) #classList has the thread ids
	#print "classList:", classList

	classL = [range(0,n),range(n,2*n),range(2*n,3*n)]
	seq_stage_id = tryDic.keys()[0] #FIXME: to be changed
	del tryDic[seq_stage_id]
	del acqDic[seq_stage_id]
	del relDic[seq_stage_id]

	end_time = maxT(relDic.values(),1)[1] #the timestamp is the second element in the tuple
	start_time = minT(tryDic.values(),1)[1]
	print "start_time:", start_time, "end_time:", end_time,"Total time:", end_time - start_time
	
	#class 1 end time
	class1EndT = start_time
	for id in (classList[0]):
		if len(relDic[id]) > 0 :
			for (lock,ts) in  relDic[id]:
				if ts > class1EndT:
					class1EndT = ts
	print "class 1 end time:", class1EndT 
	
	#print "Global analysis:"
	#res = multi_analyze(tryDic, acqDic, relDic, namesD, classL)
	#measured = res[4]
	#analyzed = mva_multiclass(*res[:4])[0]
	#print "measured class 1:", measured[:,0:1:1]
	#print "analyzed class 1:",  analyzed[1::2,0:1:1]
	#print "measured class 3:", measured[:,2::]
	#print "analyzed class 3:",  analyzed[1::2,2::]

	l = filter(lambda x: x!= [], relDic.values())
	startRel = minT(l,1)[1]
	endRel = maxT(l,1)[1]

	print "Doing analysis before class 1 threads end"
	slot_size = slotted_analyzed(start_time,class1EndT,10,0,tryDic,acqDic,relDic,namesD,n,0)
	print "Doing analysis after class 1 threads end"
	slotted_analyzed(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n,1)
	if if_print:
		print "throu1:"
		for i in throu1:
			print i
		print "throu3:"
		for i in throu3:
			print i
		print "throughput difference:"
		for i in throu_diff:
			print i
	print "n_item1"
	for i in n_item1:
		print i
	print "n_item3"
	for i in n_item3:
		print i
	print "item_diff"
	for i in item_diff:
		print i

get_derivative_analysis(10)
