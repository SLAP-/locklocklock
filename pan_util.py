#!/usr/bin/python
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
from pylab import *
import sys
from lockgraph import *

if_print = 1

class1_idx = 0
class3_idx = 2

ts_idx = 1
lock_idx = 0

start_run = 0 #the starting time of the experiment
n_locks = 57
proportions = []

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
	return a*t1+b*t1+c

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
		print "in get_bucket_proportion"
		total_bucket_size = 0
		bucket_size = 0
		for lock in range(0, n_locks):
			average_items = get_average_n_items_in_range(start,current_start,end,acqDic,lock,classList)
			if lock == lockid:
				bucket_size = average_items
			total_bucket_size = total_bucket_size + average_items
		print "bucket_size", bucket_size, "total_bucket_size:",total_bucket_size
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
	

def get_slices(tryDic,acqDic,relDic,slot_start,slot_end):
	print "in get_slices:"
	(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,slot_start,slot_end)
	newClassL = getClassL (newTryDict.keys(),tryDic.keys(),n) #print "newClassL:", newClassL
	newnewClassL = filter(lambda x: x != [], newClassL)
	classVector = tuple(map(len,newnewClassL))
	return (newTryDict, newAcqDict,newRelDict,newnewClassL)

def slotted_analyzed(start,end,nDiv,slot_size,tryDic,acqDic,relDic,namesD,n,classList):	
	print "slotted analysis:"
	throughput = []
	s1 = [] #class 1 service time
	s3 = [] #class 3 service time
	a1 = [] #analyzed class 1
	a3 = [] #analyzed class 3
	m1 = [] #measured class 1
	m3 = [] #measured class 3
	n_items = []
	serv_l = []
	ts_l = get_ts(tryDic, acqDic, relDic, slot_size, start, end, nDiv)
	for (slot_start,slot_end) in ts_l:
		(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,slot_start,slot_end)
		newClassL = getClassL (newTryDict.keys(),tryDic.keys(),n) #print "newClassL:", newClassL
		newnewClassL = filter(lambda x: x != [], newClassL)
		classVector = tuple(map(len,newnewClassL))

		serv,num_items = get_n_items_and_serv_per_slot(slot_start,slot_end,acqDic,newTryDict,newAcqDict,newRelDict,namesD,newnewClassL,n,classList)
		n_items.append(get_average_n_items_in_range(start_run,slot_start,slot_end,acqDic,0,classList))

		res = multi_analyze(newTryDict, newAcqDict, newRelDict, namesD, newClassL)
		measured = res[4][0] #first lock waiting time
		ana = mva_multiclass(res[0],res[1],classVector,res[3])
		serv_l.append(res[1])
		analyzed = ana[0] #ana[0] is the waiting time matrix, ana[1] is the queue length matrix
		fill_m_and_a(measured,analyzed,m1,m3,a1,a3)
	return serv_l, m3,a3
		
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

def iterative_throughput_s_mva(throu1, throu3, s1, s3 ,serv,res,newClassL):
	lockid = 0
	throu1 = throu1/proportions.pop(0)
	throu3 = throu3/proportions.pop(0)
	print "in iterative_throughput_s_mva: the real class 1 service time:", serv[0][class1_idx]

	for i in range(0,5): #TOFIX: change this to convergence condition
		ana = mva_multiclass(res[0],serv,newClassL,res[3]) #do an initial analysis
		(throu1,throu3) = get_throughput_from_l(ana[2][lockid])
		oh = get_contention_overhead(throu1,throu3) #overhead
		#plug the service time back in
		s1 = serv[0][class1_idx] #res_up[1] is the service time matrix
		s1 = s1 + oh
		serv[0][class1_idx] = s1 #res_up[1] is the service time matrix
		print "service 1:", s1
		print "throughput 1:", throu1
		print "throughput 3:", throu3
		#update s3 also
	
	return (ana[0][1,class1_idx],ana[0][1,class3_idx]) #ana[0] is the waiting time, ana[2] is the throughput

def predict_another_case(tryDic,acqDic,relDic,start_time,end_time,classL,newClassL,n_items_l,namesD,a1,b1,a3,b3,serv,if_up):
	predicted3 = []
	predicted1 = []
	predicted3_with_contention = []
	predicted1_with_contention = []
	print "in predict_another_case,classL:", classL, " serv length:",len(serv)
	(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,start_time,end_time)
	res = multi_analyze(newTryDict, newAcqDict, newRelDict, namesD, classL)
	if if_up:
		class3_idx = 2
	else:
		class3_idx = 0
	j = 0
	for n in n_items_l: #n is the number of items to be done
		#res[1][1][class1_idx] = 1.0/450000.7 #res_up is the service time matrix
		if j < (len(serv)):
			#for i in range(0,len(res[1])/2):
			#serv[j][i*2+1][class3_idx] = 1.0/(a*n+b) #i*2+1 are the indices of all the locks in the queueing network
			serv[j][1][class3_idx] = 1.0/(a3*n+b3) #i*2+1 are the indices of all the locks in the queueing network
			if if_up:
					#serv[j][i*2+1][class1_idx] = 1.0/(a*n+b) #res_up[1] is the service time matrix
				serv[j][1][class1_idx] = 1.0/(a1*n+b1) #res_up[1] is the service time matrix
			ana = mva_multiclass(res[0],serv[j],newClassL,res[3]) #the second parameter is the service rate, not the service time

			#ana[2][0] is the throughput
			(throu1,throu3) = get_throughput_from_l(ana[2][0])
			print "throughput:", ana[2][0] #throughput of lock 0
			if if_up:
				(w1,w3) = iterative_throughput_s_mva(throu1, throu3, serv[j][1][class1_idx], serv[j][1][class3_idx] ,serv[j],res,newClassL)
				predicted1_with_contention.append(w1)
				predicted3_with_contention.append(w3)
			predicted1.append(ana[0][1,class1_idx]) #ana[0] is the waiting time, ana[2] is the throughput
			predicted3.append(ana[0][1,class3_idx]) #ana[0] is the waiting time, ana[2] is the throughput
			#predicted3.append(w3) #ana[0] is the waiting time, ana[2] is the throughput
			j = j + 1
	#print_list_vertical(predicted3,"predicted analyzed class 3")
	return predicted1,predicted3, predicted1_with_contention,predicted3_with_contention

def pre_process(filename,n):
	print "Parsing input file:",filename
	(tryDic, acqDic, relDic, namesD) = parseDicTuple (filename)
	map(lambda x:x[0],tryDic.values())

	#classList = []
	classL = [range(0,n),range(n,2*n),range(2*n,3*n)]
	seq_stage_id = tryDic.keys()[0] #Now we assume the first thread is the sequential thread #FIXME: to be changed
	del tryDic[seq_stage_id]
	del acqDic[seq_stage_id]
	del relDic[seq_stage_id]
	threads_l = tryDic.keys()
	classList = [threads_l[0:n],threads_l[n:2*n],threads_l[2*n:3*n]]
	#print "In pre_process, classList:", classList

	end_time = maxT(relDic.values(),1)[ts_idx] #the timestamp is the second element in the tuple
	start_time = minT(tryDic.values(),1)[ts_idx]
	start_run = start_time
	#print "start_time:", start_time, "end_time:", end_time,"Total time:", end_time - start_time
	class1EndT =  get_class_end_time(start_time,classList, class1_idx, relDic) #0 is the class index
	#print "class1EndT:",class1EndT

	l = filter(lambda x: x!= [], relDic.values())
	startRel = minT(l,1)[1]
	endRel = maxT(l,1)[1]
	return tryDic,acqDic,relDic,namesD,start_time,end_time,class1EndT,classList,classL

def check_fitting(a_up,b_up,a_down,b_down,class_idx):
	n = 10
	nDiv = 10
	file_target = './dedup_run_on_halvan/dedup_10th_32c_random_input/dedup.native.10th.halvan.random'
	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(file_target,n)
	slot_size = get_slot_size(start_time,class1EndT,nDiv)

	lockid = 0
	serv_up_10, items_up_10 = get_n_items_and_serv(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n,classList,lockid)
	serv_down_10,items_down_10 =  get_n_items_and_serv(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n,classList,lockid)

	#print_list_vertical(serv_up_10,"serv_up_10")
	#print "uphill service time calculated with (a,b) of (2,2,2)"
	print "#items","\t","Real service time:", "\t","Curve fitting service time with parameters of 222"
	for index,i in enumerate(items_up_10):
		print i,"\t",serv_up_10[class_idx][index],"\t",i*a_up + b_up
	for index,i in enumerate(items_down_10):
		print i,"\t",serv_down_10[class_idx][index],"\t",i*a_down + b_down

def get_polyfit(data_file):
	a_up = a_down = b_up = b_down = 0
	n = 2
	nDiv = 10
	class1_idx = 0 
	class3_idx = 1
	lockid = 0
	print "analyzing 222"
	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(data_file,n)
	slot_size = get_slot_size(start_time,class1EndT,nDiv)
	serv_up, items_up = get_n_items_and_serv(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n,classList,lockid)

	(a1_up,b1_up) = polyfit(items_up,serv_up[class1_idx],1)
	(a3_up,b3_up) = polyfit(items_up,serv_up[class3_idx],1)
	print "result of uphill polyfit: a1_up:", a1_up, "b1_up:",b1_up
	print "result of uphill polyfit: a3_up:", a3_up, "b3_up:",b3_up

	#print "length of items_down:", len(items_down), "length of serv_down:",len(serv_down)
	serv_down,items_down =  get_n_items_and_serv(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n,classList,lockid)
	(a3_down,b3_down) = polyfit(items_down,serv_down[class3_idx],1)
	print "result of downhill polyfit:", a3_down, "b3_down:",b3_down
	#check_fitting(a_up,b_up,a_down,b_down,class_idx)
	#print "calculating the proportion of lock0 hashtable uphill:", items_up[lockid]/sum(items_up) 
	return (a1_up,b1_up,a3_up,b3_up,a3_down,b3_down)

def get_derivative_analysis():
	if_predict = int(sys.argv[2])
	file_base = './dedup_run_on_halvan/halvan.dedup.random.2th'
	file_target = './dedup_run_on_halvan/dedup_10th_32c_random_input/dedup.native.10th.halvan.random'
	n_base = 2 #num_t_base
	n_target = 10 
	nDiv = 10
	class_id = 3
	a1_down = b1_down = 0

	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(file_target,n_target)
	slot_size = get_slot_size(start_time,class1EndT,nDiv)

	serv_l_up,m3_up,a3_up = slotted_analyzed(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n_target,classList) #devide the part before class 1 threads finish into 2 parts
	serv_l_down,m3_down,a3_down = slotted_analyzed(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n_target,classList)
	print_list_vertical(m3_up + m3_down ,"measured3")
	print_list_vertical(a3_up + a3_down ,"analyzed3")

	#**************************************************************#
	#start prediction, parse the base_file#
	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(file_base,n_base)
	#get the parameters for the function (service time 3) = a * (#items in the hashtable) + b
	#here we assume a linear relationship between the number of items and the class 3 threads (compress threads) servicet item
	if  if_predict:
		(a1_up,b1_up,a3_up,b3_up,a3_down,b3_down) = get_polyfit(sys.argv[1])

		uphill_classL = (10,10,10)
		n_items_l = [216, 654, 1068, 1495, 1981, 2499, 3063, 3667, 4243, 4728]
		print "predicting uphill:"
		p1_up,p3_up,p1_w_cont_up, p3_w_cont_up = predict_another_case(tryDic,acqDic,relDic,start_time,class1EndT,classL,uphill_classL,n_items_l,namesD,a1_up,b1_up,a3_up,b3_up,serv_l_up,1)

		#print_list_vertical(predicted3_up ,"predicted3_up")
		n_items_l = [4743, 4318, 3882, 3439, 2983, 2549, 2131, 1692, 1237, 781, 275]
		print "predicting downhill:"
		downhill_classL = (10,)
		classL = [[0,1]]

		print "serv_l_down size:", len(serv_l_down)
		p1_down,p3_down,p1_w_cont_down,p3_w_cont_down =  predict_another_case(tryDic,acqDic,relDic,class1EndT,end_time,classL,downhill_classL,n_items_l,namesD,a1_down,b1_down,a3_down,b3_down,serv_l_down,0)

	print "proportions:", proportions
	#for index,i in enumerate(n_items_l):
		#print i,"\t",m3_up[index],a3_up[index],p3_up[index],p3_w_cont_up[index]
	#print len(m3_up),len(a3_up), len(p3_up),len(p3_w_cond_up)
	print "m3_up:"
	print m3_up
	print "a3_up:"
	print a3_up
	print "p3_up:"
	print p3_up
	print "p3_w_cont__up:"
	print p3_w_cont_up
	#print_list_vertical(predicted3_up + predicted3_down ,"predicted3")

get_derivative_analysis()
#analyze_measured()
#get_polyfit(sys.argv[1])
