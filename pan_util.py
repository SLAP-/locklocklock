#!/usr/bin/python
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from pylab import *
import sys
from lockgraph import *

if_print = 1

n_class = 3 #there are three classes or stages in the pipeline in our analysis
class1_idx = 0
class2_idx = 1
class3_idx = 2

ts_idx = 1
lock_idx = 0

start_run = 0 #the starting time of the experiment
n_locks = 57
proportions = []

polyfit_degree = 2

#prints out the elements of a list with a newline break after each one of them
def print_list_vertical(l,name):
	print "printing list", name	
	for i in l:
		print i

#calculate the contention overhead in the service time given t1 and t3
def get_contention_overhead(t1,t3):
	a = 8.63685e+09
	b = -1.27817e+10
	c = 56357
	return a*t1*t1+b*t1+c
def get_contention_overhead_asymmetric(t1,t3):
	a = -2.88363e+15 
	b = 8.11897e+15 
	c = -4.18053e+15
	d = 1.06434e+11 
	e = -4.45469e+10
	f = -476638    
	return a*t1*t1 + b*t3*t3 + c*t1*t3 + d*t1 + e*t3 + f

#returns new dictionaries with only the entries within the timestamp range (start_ts, end_ts)
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
			if t[1] > start_ts and t[1] <= end_ts: #the tuple index starts from 0, if the timestamp falls in the range
				newRelDict[k].append(t)
				newTryDict[k].append(tryDict[k][i])
				newAcqDict[k].append(acqDict[k][i])
			i = i + 1	
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

def get_class_list(new_threads,original_threads,n):
	newClassL = getClassL (new_threads,original_threads,n) #print "newClassL:", newClassL
	newnewClassL = filter(lambda x: x != [], newClassL)
	classVector = tuple(map(len,newnewClassL))
	return (newClassL, newnewClassL, classVector)

def count_entries(tryDict, classID, lockID, start_ts, end_ts,classList):
	count = 0
	for k,v in tryDict.iteritems(): #v is a list of the tuple (lockid,timstamp)
		if k in classList[classID]:
			for x in v:
				if x[0] == lockID and x[1] > start_ts and x[1] <= end_ts:
					count = count + 1
	return count 
	
def fill_m_and_a(measured,analyzed,m1,m3,a1,a3):
	if len(measured) == 3:
		m1.append(measured[0])
		m3.append(measured[2])
		a1.append(analyzed[1][0])
		a3.append(analyzed[1][2])
	elif len(measured) == 2:
		m1.append(measured[0])
		m3.append(measured[1])
		a1.append(analyzed[1][0])
		a3.append(analyzed[1][1])
	elif len(measured) == 1:
		m3.append(measured[0])
		a3.append(analyzed[1][0])

def get_average_n_items_in_range(start,current_start,end,acqDic,lockid,classList):
		prev_put = count_entries(acqDic,class1_idx,lockid,start,current_start,classList)
		prev_got = count_entries(acqDic,class3_idx,lockid,start,current_start,classList)
		n_begin = prev_put - prev_got
		n1 = count_entries(acqDic,class1_idx,lockid,current_start,end,classList)
		n3 = count_entries(acqDic,class3_idx,lockid,current_start,end,classList)
		return  n_begin + (n1 - n3)/2

def get_n_items_and_serv_per_slot(slot_start,slot_end,acqDic,newTryDict,newAcqDict,newRelDict,namesD,newnewClassL,n,classList):
	serv = multi_serviceTime(newTryDict,newAcqDict,newRelDict,namesD,newnewClassL,0) #0 is the overhead
	num_items = np.empty(n_locks,list)
	for i in range(0,n_locks):
		num_items[i] = []
	for i in range(0,n_locks): #for all the locks
		num_items[i].append(get_average_n_items_in_range(start_run,slot_start,slot_end,acqDic,i,classList))
	return serv,num_items

def get_slot_size(start,end,nDiv):
	assert nDiv > 0
	return (end-start)/nDiv

def get_num_div(start,end,slot_size):
	assert slot_size > 0
	return (end-start)/slot_size

def get_bucket_proportion(start,current_start,end,acqDic,lockid,classList):
		#print "in get_bucket_proportion, classList", classList
		total_bucket_size = 0
		bucket_size = 0
		for lock in range(0, n_locks):
			average_items = get_average_n_items_in_range(start,current_start,end,acqDic,lock,classList)
			#print "lock:",lock, "aveage_items:",average_items
			if lock == lockid:
				bucket_size = average_items
			total_bucket_size = total_bucket_size + average_items
		#print "bucket_size", bucket_size, "total_bucket_size:",total_bucket_size
		return (float(bucket_size))/total_bucket_size
		
def get_n_items_and_serv(start,end,nDiv,slot_size,tryDic,acqDic,relDic,namesD,n,classList,lockid):
	service1 = [] #class 1 service time
	service3 = []
	n_items = np.empty(n_locks,list)
	for i in range(0,len(n_items)):
		n_items[i] = []
	ts_l = get_ts(tryDic, acqDic, relDic, slot_size, start, end, nDiv)
	for (slot_start,slot_end) in ts_l:
		(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,slot_start,slot_end)
		newClassL = getClassL (newTryDict.keys(),tryDic.keys(),n) 
		newnewClassL = filter(lambda x: x != [], newClassL)
		serv,num_items = get_n_items_and_serv_per_slot(slot_start,slot_end,acqDic,newTryDict,newAcqDict,newRelDict,namesD,newnewClassL,n,classList)
		for i in range(0,len(n_items)): #for all the locks
			n_items[i].append(get_average_n_items_in_range(start_run,slot_start,slot_end,acqDic,i,classList))
		if len(serv) == 3:		
			service1.append(serv[0][0])
			service3.append(serv[2][0])
		elif len(serv) == 2:
			service1.append(serv[0][0])
			service3.append(serv[1][0])
		elif len(serv) == 1:
			service3.append(serv[0][0])
		prop = get_bucket_proportion(start_run,slot_start,slot_end,acqDic,0,classList)
		proportions.append(prop)
	return [service1, service3], n_items[lockid]

def get_ts(tryDic,acqDic,relDic,slot_size, start, end, nDiv):
	ts_l = []
	if slot_size == 0 and nDiv > 0:
		slot_size = get_slot_size(start,end,nDiv)
	elif nDiv == 0:
		nDiv = get_num_div(start,end,slot_size)

	for i in range(nDiv):
		slot_start = start + i*slot_size 
		if i == nDiv - 1:
			slot_end = end
		else:
			slot_end = slot_start + slot_size
		ts_l.append((slot_start,slot_end))
	return ts_l
	
def slotted_analyzed(start,end,nDiv,slot_size,tryDic,acqDic,relDic,namesD,n,classList):	
	#print "slotted analysis: classList:", classList
	throughput = []
	s1 = [] #class 1 service time
	s3 = [] #class 3 service time
	a1 = [] #analyzed class 1
	a3 = [] #analyzed class 3
	m1 = [] #measured class 1
	m3 = [] #measured class 3
	lockid = 0
	n_items = []
	serv_l = []
	ts_l = get_ts(tryDic, acqDic, relDic, slot_size, start, end, nDiv)
	for (slot_start,slot_end) in ts_l:
		(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,slot_start,slot_end)
		(newClassL,newnewClassL, classVector) = get_class_list(newTryDict.keys(),tryDic.keys(),n)

		serv,num_items = get_n_items_and_serv_per_slot(slot_start,slot_end,acqDic,newTryDict,newAcqDict,newRelDict,namesD,newnewClassL,n,classList)
		n_items.append(num_items[lockid])
		res = multi_analyze(newTryDict, newAcqDict, newRelDict, namesD, newClassL)
		measured = res[4][0] #first lock waiting time
		ana = mva_multiclass(res[0],res[1],classVector,res[3])
		serv_l.append(res[1])
		analyzed = ana[0] #ana[0] is the waiting time matrix, ana[1] is the queue length matrix
		fill_m_and_a(measured,analyzed,m1,m3,a1,a3)
	#print "in slotted analysis, printing a1:", a1
	#print "in slotted analyzed: num_items:", num_items
	return serv_l, n_items,m3,a3,m1,a1
		
def get_class_end_time(start_time,classList,classID, relDic):
	end_time = start_time
	for id in classList[classID]:
		if len(relDic[id]) > 0:
			for(lock,ts) in relDic[id]:
				if ts > end_time:
					end_time = ts
	return end_time

#given a list of throughputs, return the tuple (class1_throughput, class3_throughput)
def get_throughput_from_l(l):
	if len(l) == 3:
		return (l[0],l[2])
	elif len(l) == 2:
		return (l[0],l[1])
	elif len(l) == 0: #class 1 has already disappeared
		return (0,l[0])
	else:
		return (0,l)

def iterative_throughput_s_mva(throu1, throu3, s1, s3 ,serv,res,newClassL): #s1 and s3 and serv is service rate
	lockid = 0
	prop = proportions.pop(0)
	i = 0
	#print "in iterative_throughput_s_mva: the proportion:", prop
	print "in iterative_throughput_s_mva: the real class 1 service rate:", serv[0][class1_idx]
	print "in iterative_throughput_s_mva: the real class 1 real service time:", 1.0/serv[0][class1_idx]
	print "in iterative_throughput_s_mva: the polyfit class 1 service time:", 1.0/s1
	
	#if if_print: #TOFIX: change this to convergence condition
	for i in range(0,5): #TOFIX: change this to convergence condition
		ana = mva_multiclass(res[0],serv,newClassL,res[3]) #do an initial analysis
		(throu1,throu3) = get_throughput_from_l(ana[2][lockid])
		throu1_stage = throu1/prop
		throu3_stage = throu3/prop
		print "throughput 1:", throu1_stage, "throughput 3:",throu3_stage
	
		oh = get_contention_overhead(throu1_stage,throu3_stage) #overhead
		oh_asy = get_contention_overhead_asymmetric(throu1_stage,throu3_stage) #overhead
		#plug the service time back in
		s1 = 1.0/serv[0][class1_idx] #res_up[1] is the service time matrix
		print "service 1 before plugging overhead:", s1
		s1 = s1 + oh_asy
		serv[0][class1_idx] = 1.0/s1 #res_up[1] is the service time matrix
		print "overhead with the asymmetric formula:", oh_asy
		print "service 1:", s1
		i = i + 1
	return (ana[0][1,class1_idx],ana[0][1,class3_idx]) #ana[0] is the waiting time, ana[2] is the throughput

def predict_another_case(tryDic,acqDic,relDic,start_time,end_time,classL,n_items_l,namesD,fit_param,serv,if_up,n_base,n_target): #serv has the service rates, not time
	print "in predict another case:"
	print "classL:", classL
	predicted3 = []
	predicted1 = []
	predicted3_with_contention = []
	predicted1_with_contention = []
	a1 = fit_param[0][0][0]
	a3 = fit_param[1][0][0]
	b1 = fit_param[0][0][1]
	b3 = fit_param[1][0][1]

	(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,start_time,end_time)
	(newClassL,newnewClassL,classVector) = get_class_list(newTryDict.keys(),tryDic.keys(),n_base)
	print "newClassL:", newClassL
	print "newnewClassL:", newnewClassL
	print "classVector:", classVector
	res = multi_analyze(newTryDict, newAcqDict, newRelDict, namesD, classL)
	if if_up:
		class3_idx = 2
	else:
		class3_idx = 0
	for n in n_items_l: #n is the number of items to be done
		serv[1][class3_idx] = 1.0/(a3*n+b3) #i*2+1 are the indices of all the locks in the queueing network
		if if_up:
			serv[1][class1_idx] = 1.0/(a1*n+b1) #res_up[1] is the service time matrix
		#ana = mva_multiclass(res[0],serv,newClassL,res[3]) #the second parameter is the service rate, not the service time
		#ana = mva_multiclass(res[0],res[1],newClassL,res[3])
		print "res[2]:",res[2]
		#ana = mva_multiclass(res[0],serv,newClassL,res[3])
		ana = mva_multiclass(res[0],serv,(10,10,10),res[3])
		#ana = mva_multiclass(res[0],serv,res[2],res[3])
		#(throu1,throu3) = get_throughput_from_l(ana[2][0])
		#print "throughput before calling iterative throughput", ana[2][0] #throughput of lock 0
		#print "number of items:", n
		#if if_up:
			#(w1,w3) = iterative_throughput_s_mva(throu1, throu3, serv[j][1][class1_idx], serv[j][1][class3_idx], serv[j], res, newClassL)
			#predicted1_with_contention.append(w1)
			#predicted3_with_contention.append(w3)
		predicted1.append(ana[0][1,class1_idx]) #ana[0] is the waiting time, ana[2] is the throughput
		predicted3.append(ana[0][1,class3_idx]) #ana[0] is the waiting time, ana[2] is the throughput
			#predicted3.append(w3) #ana[0] is the waiting time, ana[2] is the throughput
	return predicted1,predicted3, predicted1_with_contention,predicted3_with_contention

def pre_process(filename,n):
	print "Parsing input file:",filename
	(tryDic, acqDic, relDic, namesD) = parseDicTuple (filename)
	map(lambda x:x[0],tryDic.values())

	classL = [range(0,n),range(n,2*n),range(2*n,3*n)]
	seq_stage_id = tryDic.keys()[0] #Now we assume the first thread is the sequential thread #FIXME: to be changed
	class_end_time = [0 for i in range(n_class)]
	del tryDic[seq_stage_id]
	del acqDic[seq_stage_id]
	del relDic[seq_stage_id]
	threads_l = tryDic.keys()
	classList = [threads_l[0:n],threads_l[n:2*n],threads_l[2*n:3*n]]

	end_time = maxT(relDic.values(),1)[ts_idx] #the timestamp is the second element in the tuple
	start_time = minT(tryDic.values(),1)[ts_idx]
	start_run = start_time
	#print "start_time:", start_time, "end_time:", end_time,"Total time:", end_time - start_time
#	class1EndT = get_class_end_time(start_time,classList,class1_idx, relDic) #0 is the class index
#	class2EndT = get_class_end_time(start_time,classList,class2_idx, relDic)
#	class3EndT = get_class_end_time(start_time,classList,class3_idx, relDic)
	for i in range(n_class):
		class_end_time[i] = get_class_end_time(start_time,classList,i, relDic)
		print "class ", i, "end time:", class_end_time[i]
	l = filter(lambda x: x!= [], relDic.values())
	startRel = minT(l,1)[1]
	endRel = maxT(l,1)[1]
	return tryDic,acqDic,relDic,namesD,start_time,end_time,class_end_time,classList,classL

def get_polyfit(tryDic,acqDic,relDic,namesD,start_time,end_time,stage_end_time,classList,classL):
	a_up = a_down = b_up = b_down = 0
	#n = 2
	n = 4
	nDiv = 10
	class1_idx = 0 
	class3_idx = 1
	lockid = 0
	n_ts_chunks = 3
	a_idx = 0
	b_idx = 1
	serv = [ [ [] for i in range(n_ts_chunks)] for j in range(n_class)]
	fit_param = [[ [] for i in range(n_ts_chunks)] for j in range(n_class)]

	class1EndT = stage_end_time[class1_idx]
	class2EndT = stage_end_time[class2_idx]
	class3EndT = stage_end_time[class3_idx]

	slot_size = get_slot_size(start_time,class1EndT,nDiv)

	print "between the beginning and class 2 ends:"
	[serv[class1_idx][0],serv[class3_idx][0]], items = get_n_items_and_serv(start_time,class2EndT,10,0,tryDic,acqDic,relDic,namesD,n,classList,lockid)
	(a1_up,b1_up) = polyfit(items,serv[class1_idx][0],1)
	(a,b) = polyfit(items,serv[class1_idx][0],1)
	fit_param[class1_idx][0].append(a)
	fit_param[class1_idx][0].append(b)

	(a,b) = polyfit(items,serv[class3_idx][0],1)
	fit_param[class3_idx][0].append(a)
	fit_param[class3_idx][0].append(b)

	[serv[class1_idx][1],serv[class3_idx][1]], items = get_n_items_and_serv(class2EndT,class1EndT,10,0,tryDic,acqDic,relDic,namesD,n,classList,lockid)

	print "between class 2 ends and class 1 ends:"
	(a,b) = polyfit(items,serv[class1_idx][1],1)
	fit_param[class1_idx][1].append(a)
	fit_param[class1_idx][1].append(b)

	(a,b) = polyfit(items,serv[class3_idx][1],1)
	fit_param[class3_idx][1].append(a)
	fit_param[class3_idx][1].append(b)

	serv_down,items_down =  get_n_items_and_serv(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n,classList,lockid)
	(a,b) = polyfit(items_down,serv_down[class3_idx],1)
	fit_param[class3_idx][2].append(a)
	fit_param[class3_idx][2].append(b)

	print_list_vertical(items_down,"number of items down")
	print_list_vertical(serv_down[class3_idx],"service time class 3 down")
	#return (a1_up,b1_up,a3_up,b3_up,a3_down,b3_down)
	return fit_param

def get_derivative_analysis():
	if_predict = int(sys.argv[1])
	#file_base = './dedup_run_on_halvan/halvan.dedup.random.2th'
	file_base = './dedup_run_on_halvan/halvan.dedup.random.4th'
	file_target = './dedup_run_on_halvan/dedup_10th_32c_random_input/dedup.native.10th.halvan.random'
	#n_base = 2 #num_t_base
	n_base = 4
	n_target = 10 
	nDiv = 10
	n_items = []

	tryDic,acqDic,relDic,namesD,start_time,end_time, stage_end_time ,classList,classL = pre_process(file_target,n_target)
	class1EndT = stage_end_time[class1_idx]
        class2EndT = stage_end_time[class2_idx]
        class3EndT = stage_end_time[class3_idx]

	slot_size = get_slot_size(start_time,class1EndT,nDiv)
	
	s = []
	m1_total = []
	m3_total = []
	a1_total = []
	a3_total = []

	serv_l_up,num_items,m3,a3,m1,a1 = slotted_analyzed(start_time,class2EndT,1,0,tryDic,acqDic,relDic,namesD,n_target,classList) #devide the part before class 1 threads finish into 2 parts
	#print "length of serv_l_up:",len(serv_l_up)
	m1_total.append(m1)
	a1_total.append(a1)
	m3_total.append(m3)
	a3_total.append(a3)
	print "num_items: 1st ", num_items 
	n_items.append(num_items)

	serv_l_up2,num_items,m3,a3,m1,a1 = slotted_analyzed(class2EndT,class1EndT,1,0,tryDic,acqDic,relDic,namesD,n_target,classList) #devide the part before class 1 threads finish into 2 parts
	m1_total.append(m1)
	a1_total.append(a1)
	m3_total.append(m3)
	a3_total.append(a3)
	n_items.append(num_items)

	serv_l_down,num_items,m3,a3,m1,a1 = slotted_analyzed(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n_target,classList)

	m1_total.append(m1)
	a1_total.append(a1)
	m3_total.append(m3)
	a3_total.append(a3)
	n_items.append(num_items)

	print "m1:",m1_total
	print "a1:",a1_total
	print "m3:",m3_total
	print "a3:",a3_total

	#**************************************************************#
	#get the parameters for the function (service time 3) = a * (#items in the hashtable) + b
	#here we assume a linear relationship between the number of items and the class 3 threads (compress threads) servicet item
	if if_predict:
		tryDic,acqDic,relDic,namesD,start_time,end_time, stage_end_time, classList,classL = pre_process(file_base,n_base)
		#print "in if_predict:"
		#¤print "classList:", classList
		#print "classL:", classL
		class1EndT = stage_end_time[class1_idx]
        	class2EndT = stage_end_time[class2_idx]
        	class3EndT = stage_end_time[class3_idx]

		fit_param = get_polyfit(tryDic,acqDic,relDic,namesD,start_time,end_time,stage_end_time,classList,classL)
		#print "fit_param:",fit_param
		#print "num_items:", n_items
		
		#uphill_classL = (10,10,10)
		#n_items_l = [216, 654, 1068, 1495, 1981, 2499, 3063, 3667, 4243, 4728]
		#print "predicting uphill:"
		#print "serv_l_up:"
		#print serv_l_up
		p1_up,p3_up,p1_w_cont_up, p3_w_cont_up = predict_another_case(tryDic,acqDic,relDic,start_time,class2EndT,classL,num_items[0],namesD,fit_param,serv_l_up[0],1,n_base,n_target)
		#n_items_l = [4743, 4318, 3882, 3439, 2983, 2549, 2131, 1692, 1237, 781, 275]
		print "p1:", p1_up
		#print "predicting downhill:"
		#downhill_classL = (10,)
		#classL = [[0,1]]
get_derivative_analysis()
