#!/usr/bin/python
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
import sys
from lockgraph import *

if_print = 1
classList = []
n_item1 = []
n_item3 = []
item_diff = []

def print_list_vertical(l,name):
	print "printing list", name	
	for i in l:
		print i

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
			if t[1] > start_ts and t[1] <= end_ts: #the tuple index starts from 0, if the timestamp falls in the rage
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
	
def slotted_analyzed(start,end,nDiv,slot_size,tryDic,acqDic,relDic,namesD,n):
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
		slot_start = start + i*slot_size 
		if i == nDiv - 1:
			slot_end = end
		else:
			slot_end = slot_start + slot_size

		#print "slot starts: ",slot_start, " ends: ",slot_end, " slot size:", slot_end - slot_start
		(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,slot_start,slot_end)
		newClassL = getClassL (newTryDict.keys(),tryDic.keys(),n)
		#print "newClassL:", newClassL

		res = multi_analyze(newTryDict, newAcqDict, newRelDict, namesD, newClassL)
		measured = res[4][0] #first lock waiting time
		#print "measured:", measured
		newnewClassL = filter(lambda x: x != [], newClassL)
		#print "newClassL:", newClassL

		serv = multi_serviceTime(newTryDict,newAcqDict,newRelDict,namesD,newnewClassL,0)
		#print "newnewClassL:", newnewClassL
		classVector = tuple(map(len,newnewClassL))
		#analyzed = mva_multiclass(*res[:4])[0]
		print "classVector:", classVector
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
	print_list_vertical(a3,"analyzed class 3")

	if len(s3) == len(item_diff):
		return slot_size, s3, item_diff
	else:
		return slot_size, s3, item_diff[len(item_diff) - len(s3):]
		

def get_class_end_time(start_time,classList,classID, relDic):
	end_time = start_time
	for id in classList[classID]:
		if len(relDic[id]) > 0:
			for(lock,ts) in relDic[id]:
				if ts > end_time:
					end_time = ts
	return end_time
def global_mva_analysis(tryDic, acqDic,relDic,namesD,classL):
	print "Global analysis:"
	res = multi_analyze(tryDic, acqDic, relDic, namesD, classL)
	measured = res[4]
	analyzed = mva_multiclass(*res[:4])[0]
	return res

def get_derivative_analysis(n):
	class1_idx = 0
	class3_idx = 2
	ts_idx = 1
	lock_idx = 0
	print "Parsing input file:",sys.argv[1]
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

	end_time = maxT(relDic.values(),1)[ts_idx] #the timestamp is the second element in the tuple
	start_time = minT(tryDic.values(),1)[ts_idx]
	#print "start_time:", start_time, "end_time:", end_time,"Total time:", end_time - start_time
	class1EndT =  get_class_end_time(start_time,classList, class1_idx, relDic) #0 is the class index

	l = filter(lambda x: x!= [], relDic.values())
	startRel = minT(l,1)[1]
	endRel = maxT(l,1)[1]

	print "Doing analysis before class 1 threads end"
	#slot_size = slotted_analyzed(start_time,class1EndT,10,0,tryDic,acqDic,relDic,namesD,n,0)
	slot_size,service3_up,n_token_up = slotted_analyzed(start_time,class1EndT,10,0,tryDic,acqDic,relDic,namesD,n)
	(a_up,b_up) = polyfit(n_token_up,service3_up,1)
	print "result of polyfit: a_up:", a_up, "b_up:",b_up
	print "Doing analysis after class 1 threads end"
	sii,service3_down,n_token_down = slotted_analyzed(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n)
	(a_down,b_down) = polyfit(n_token_down,service3_down,1)
	print "result of downhill polyfit of the (amount of work, class 3 service time): a_down:", a_down, "b_down:",b_down
	predicted3 = []
	uphill_classL = (10,10,10)
	(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,start_time,class1EndT)
	res_up =  global_mva_analysis(newTryDict, newAcqDict,newRelDict,namesD,classL) #res_up returns tuple ()
	for n in (393,798,1152,1578,2032,2517,3050,3604,4093,4477): #n is the number of items to be done
	#Do the uphill analysis for the different number of items left in the queue
		#print "at the end classL:", classL
		#print "service rate of lock 0 for 222 case:", res_up[1][lock_idx]
		#print "newly calculated service time for lock 0:", (a_up*n+b_up)
		#if we recalculate the service time for all the locks with the linear equation
		print "the size of res_up[1]:", len(res_up[1])
		for i in range(0,len(res_up[1])/2):
			#res_up[1][lock_idx][class3_idx] = 1.0/(a_up*n+b_up) #res_up is the service time matrix
			res_up[1][i*2+1][class3_idx] = 1.0/(a_up*n+b_up) #res_up is the service time matrix
		#res_up[1][lock_idx][class1_idx] = 1.0/(a_up*n+b_up)
		ana = mva_multiclass(res_up[0],res_up[1],uphill_classL,res_up[3]) #the second parameter is the service rate, not the service time
		print "analyzed class 3:",  ana[0][1,class3_idx] #1 is the place of the first lock	#ana[0] is the waiting time matrix, [1::2]: gets every other row
		predicted3.append(ana[0][1,class3_idx])
	print_list_vertical(predicted3,"predicted class 3")
get_derivative_analysis(2)
