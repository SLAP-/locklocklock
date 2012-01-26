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

def print_list_vertical(l,name):
	print "printing list", name	
	for i in l:
		print i

def get_n_threads_from_dict(d):
	return d.keys()-1

def plot_data(x,y,tit,x_label,y_label,color):
	title(tit)
	xlabel(x_label)
	ylabel(y_label)
	if color == 0:
		plot(x,y,'r.')
	else:
		plot(x,y,'y.')
	show()

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

def get_n_items_and_serv(start,end,nDiv,slot_size,tryDic,acqDic,relDic,namesD,n,classList,class_id):
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
	if class_id == 1:
		return service1, n_items[0]
	else:
		return service3, n_items[0]

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
	a1 = []
	a3 = []
	m1 = []
	m3 = []
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

	print_list_vertical(m1,"measured class 1")
	print_list_vertical(m3,"measured class 3")
	print_list_vertical(a1,"analyzed class 1")
	print_list_vertical(a3,"analyzed class 3")
	print_list_vertical(n_items,"number of items")
	return serv_l
		
def get_class_end_time(start_time,classList,classID, relDic):
	end_time = start_time
	for id in classList[classID]:
		if len(relDic[id]) > 0:
			for(lock,ts) in relDic[id]:
				if ts > end_time:
					end_time = ts
	return end_time

def predict_another_case(tryDic,acqDic,relDic,start_time,end_time,classL,newClassL,n_items_l,namesD,a,b,serv,if_up):
	predicted3 = []
	print "in predict_another_case,classL:", classL, " serv length:",len(serv)
	(newTryDict,newAcqDict,newRelDict) = keep_ts_in_range(tryDic,acqDic,relDic,start_time,end_time)
	res = multi_analyze(newTryDict, newAcqDict, newRelDict, namesD, classL)
	if if_up:
		class3_idx = 2
	else:
		class3_idx = 0
	j = 0
	for n in n_items_l: #n is the number of items to be done
		#if j > 5:
			#classL = [[0,1],[],[2,3]]
			#newClassL = (10,0,10)
		#print "serv[j] length:",len(serv[j])
		#if we recalculate the service time for all the locks with the linear equation
		#for i in range(0,len(res[1])/2):
			#for k in range(0, len(res[1][i])):
			#print "res[1][i] length:", len(res[1][i]), "serv[j][i] length:",len(serv[j][i])
			#res[1][i][class1_idx] = 1.0/serv[j][i][class1_idx] #res[1] is the service time matrix
			#res[1][i][class3_idx] = 1.0/serv[j][i][class3_idx] #res[1] is the service time matrix
		#	res[1][i*2+1][class3_idx] = 1.0/((a*n*sum(n_items[i])/sum(n_items[0]))+b) #res_up[1] is the service time matrix
		#res[1][1][class1_idx] = 1.0/450000.7 #res_up is the service time matrix
		if j < (len(serv)):
			print "in predict_another_case,j:", j
			for i in range(0,len(res[1])/2):
				serv[j][i*2+1][class3_idx] = 1.0/(a*n+b) #res_up[1] is the service time matrix
			ana = mva_multiclass(res[0],serv[j],newClassL,res[3]) #the second parameter is the service rate, not the service time
			
			print "analyzed class 3:",  ana[0][1,class3_idx] #1 is the place of the first lock	#ana[0] is the waiting time matrix, [1::2]: gets every other row
			predicted3.append(ana[0][1,class3_idx])
		#throu.append(ana[2][0])
			j = j + 1
	#return throu #ana[2] is the throughput, index 0 shows the first lock
	print_list_vertical(predicted3,"predicted analyzed class 3")

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
	print "In pre_process, classList:", classList

	end_time = maxT(relDic.values(),1)[ts_idx] #the timestamp is the second element in the tuple
	start_time = minT(tryDic.values(),1)[ts_idx]
	start_run = start_time
	print "start_time:", start_time, "end_time:", end_time,"Total time:", end_time - start_time
	class1EndT =  get_class_end_time(start_time,classList, class1_idx, relDic) #0 is the class index
	print "class1EndT:",class1EndT

	l = filter(lambda x: x!= [], relDic.values())
	startRel = minT(l,1)[1]
	endRel = maxT(l,1)[1]
	return tryDic,acqDic,relDic,namesD,start_time,end_time,class1EndT,classList,classL

def list_minus(l1,l2):
	l = []
	for i in l1:
		if not i in l2:
			l.append(i)
	return l		

def get_polyfit(data_file):
	n = 2
	nDiv = 10
	class_id = 3
	print "analyzing 222"
	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(data_file,n)
	slot_size = get_slot_size(start_time,class1EndT,nDiv)
	serv_up, items_up = get_n_items_and_serv(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n,classList,class_id)

	(a_up,b_up) = polyfit(items_up,serv_up,1)
	print "result of uphill polyfit: a_up:", a_up, "b_up:",b_up
	serv_down,items_down =  get_n_items_and_serv(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n,classList,class_id)
	print "length of items_down:", len(items_down), "length of serv_down:",len(serv_down)
	(a_down,b_down) = polyfit(items_down,serv_down,1)
	print "result of downhill polyfit:", a_down, "b_down:",b_down

	file_target = './dedup_run_on_halvan/dedup_10th_32c_random_input/dedup.native.10th.halvan.random'
	n = 10
	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(file_target,n)
	slot_size = get_slot_size(start_time,class1EndT,nDiv)

	print "analyzing 101010"
	serv_up_10, items_up_10 = get_n_items_and_serv(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n,classList,class_id)

	serv_down_10,items_down_10 =  get_n_items_and_serv(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n,classList,class_id)
	print "in get_polyfit:"
	print_list_vertical(serv_up_10,"serv_up_10")
	print "uphill service time calculated with (a,b) of (2,2,2)"
	for i in items_up_10:
		print i*a_up + b_up
	
	print "downhill real service time:"
	print_list_vertical(serv_down_10,"serv_down_10")
	print "downhill service time calculated with (a,b) of (2,2,2)"
	for i in items_down_10:
		print i*a_down + b_down 
	return (a_up,b_up,a_down,b_down)

def get_derivative_analysis():
	if_predict = int(sys.argv[2])
	file_base = './dedup_run_on_halvan/halvan.dedup.random.2th'
	file_target = './dedup_run_on_halvan/dedup_10th_32c_random_input/dedup.native.10th.halvan.random'
	n_base = 2 #num_t_base
	n_target = 10 
	nDiv = 10
	class_id = 3

	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(file_target,n_target)
	slot_size = get_slot_size(start_time,class1EndT,nDiv)

	serv_l_up = slotted_analyzed(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n_target,classList) #devide the part before class 1 threads finish into 2 parts
	serv_l_down = slotted_analyzed(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n_target,classList)

	#**************************************************************#
	#start prediction, parse the base_file#
	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(file_base,n_base)
	#get the parameters for the function (service time 3) = a * (#items in the hashtable) + b
	#here we assume a linear relationship between the number of items and the class 3 threads (compress threads) servicet item
	if  if_predict:
		(a_up,b_up,a_down,b_down) = get_polyfit(sys.argv[1])

		uphill_classL = (10,10,10)
		n_items_l = [216, 654, 1068, 1495, 1981, 2499, 3063, 3667, 4243, 4728]
		print "predicting uphill:"
		throughput_up = predict_another_case(tryDic,acqDic,relDic,start_time,class1EndT,classL,uphill_classL,n_items_l,namesD,a_up,b_up,serv_l_up,1)

		n_items_l = [4743, 4318, 3882, 3439, 2983, 2549, 2131, 1692, 1237, 781, 275]
		print "predicting downhill:"
		downhill_classL = (10,)
		classL = [[0,1]]

		print "serv_l_down size:", len(serv_l_down)
		throughput_down = predict_another_case(tryDic,acqDic,relDic,class1EndT,end_time,classL,downhill_classL,n_items_l,namesD,a_down,b_down,serv_l_down,0)

#def analyze_measured():
	#initlize the number of items left and put and got arrays
#	data_file = sys.argv[1]
#	n = (int)sys.argv[2]
#	nDiv = 10
#	tryDic,acqDic,relDic,namesD,start_time,end_time, class1EndT,classList,classL = pre_process(data_file,n)
#	slot_size = get_slot_size(start_time,class1EndT,nDiv)

#	serv_l_up = slotted_analyzed(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n_target,classList) #devide the part before class 1 threads finish into 2 parts
#	print "calling get_n_items_and_serv with up hill parameters:"
#	serv_up, items_up = get_n_items_and_serv(start_time,class1EndT,2,0,tryDic,acqDic,relDic,namesD,n_target,classList)
#	s_u = serv_up[:]
#	i_u = items_up[:]
#	serv_l_down = slotted_analyzed(class1EndT,end_time,0,slot_size,tryDic,acqDic,relDic,namesD,n_target,classList)




get_derivative_analysis()
#analyze_measured()
#get_polyfit(sys.argv[1])
