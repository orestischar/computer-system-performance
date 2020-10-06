import math
import numpy as np
import sys
import random
import bisect

cpu = {"meanService_s":24, "meanService_v":21, "oldClock":0.0, "sumTimeLength":0.0, "sumBusyTime":0.0, "numberCompletions":0, "numberCompletions_s":0, "numberCompletions_v":0, "numberServers":1, "bt":0.0, "tl":0.0, "nc":0.0, "nc_v":0.0, "nc_s":0.0, "ncsq":0, "ncsq_s":0, "tlsq":0.0, "tlxcl":0.0, "tlxnc":0.0, "varnc":0.0, "covartlnc":0, "dql":0, "qt":0, "ql":0, "dqt":0, "dqt_s":10, "vartl":0}
disk = {"meanService_s":25, "meanService_v":30, "oldClock":0.0, "sumTimeLength":0.0, "sumBusyTime":0.0, "numberCompletions":0, "numberCompletions_s":0, "numberCompletions_v":0, "numberServers":1, "bt":0.0, "tl":0.0, "nc":0.0, "nc_v":0.0, "nc_s":0.0, "ncsq":0, "ncsq_s":0, "tlsq":0.0, "tlxcl":0.0, "tlxnc":0.0, "varnc":0.0, "covartlnc":0, "dql":0, "qt":0, "ql":0, "dqt":0, "dqt_s":10, "vartl":0}
out = {"meanService_s":436, "meanService_v":321, "oldClock":0.0, "sumTimeLength":0.0, "sumBusyTime":0.0, "numberCompletions":0, "numberCompletions_s":0, "numberCompletions_v":0, "numberServers":1, "bt":0.0, "tl":0.0, "nc":0.0, "nc_v":0.0, "nc_s":0.0, "ncsq":0, "ncsq_s":0, "tlsq":0.0, "tlxcl":0.0, "tlxnc":0.0, "varnc":0.0, "covartlnc":0, "dql":0, "qt":0, "ql":0, "dqt":0, "dqt_s":10, "vartl":0}

happening = [math.inf, math.inf, math.inf, math.inf, math.inf]
#[cpu, disk, balkers, outers, new_arrival]


cQ=[]
dQ=[]
oQ=[]
bQ=[]

meandiskvisits_s=17
meandiskvisits_v=9


def updateQueue(station, tup):
	global clock1, happening, cQ, dQ, oQ, bQ, subscribers_in_system, visitors_in_system, N_s, N_v

	if(station==0): #cpu

		cpu["sumTimeLength"] = cpu["sumTimeLength"] + (clock1 - cpu["oldClock"])*len(cQ)
		cpu["sumBusyTime"] = cpu["sumBusyTime"] + (clock1 - cpu["oldClock"])*min(len(cQ),cpu["numberServers"])
		cpu["oldClock"] = clock1
		
		bisect.insort(cQ, tup)
		happening[0] = clock1+cQ[0][0]


	elif(station==1): #disk

		disk["sumTimeLength"] = disk["sumTimeLength"] + (clock1 - disk["oldClock"])*len(dQ)
		disk["sumBusyTime"] = disk["sumBusyTime"] + (clock1 - disk["oldClock"])*min(len(dQ),disk["numberServers"])
		disk["oldClock"] = clock1

		dQ.append(tup)
		happening[1] = clock1+dQ[0][0]

	elif(station==2):
		if(tup[4]=='V'):
			tup_b=[ tup[3], tup[1] ]
			bisect.insort(bQ, tup_b)
			happening[2] = clock1+bQ[0][0]

	elif(station==3): #out
		out["sumTimeLength"] = out["sumTimeLength"] + (clock1 - out["oldClock"])*len(oQ)
		out["sumBusyTime"] = out["sumBusyTime"] + (clock1 - out["oldClock"])*min(len(oQ),out["numberServers"])
		out["oldClock"] = clock1

		oQ.append(tup)
		happening[3] = clock1+oQ[0][0]

def minf(happening):
	minn=happening[0]
	gurna=0
	for i,j in zip(happening, range(len(happening))):		
		if(i<minn):
			minn=i
			gurna=j
	return gurna


def ReductDatTime(time):
	global happening, cQ, dQ, oQ, bQ

	for i in range(len(cQ)):
		cQ[i][0]-=time
	for i in range(len(bQ)):
		bQ[i][0]-=time
	if dQ:
		dQ[0][0]-=time
	if oQ:
		oQ[0][0]-=time

def pop_it(index):
	global happening, cQ, dQ, oQ, bQ, clock1

	flag=True

	for i in range(len(cQ)):
		if cQ[i][1]==index:
			cQ.pop(i)
			if i==0:
				if cQ:
					happening[0]=clock1+cQ[0][0]
				else: 
					happening[0]=math.inf
			flag=False
			break
	if (flag):
		for i in range(len(dQ)):
			if dQ[i][1]==index:
				dQ.pop(i)
				if i==0:
					if dQ:
						happening[1]=clock1+dQ[0][0]
					else:
						happening[1]=math.inf
				flag=False
				break
	if (flag):
		for i in range(len(oQ)):
			if oQ[i][1]==index:
				oQ.pop(i)
				if i==0:
					if oQ:
						happening[3]=clock1+oQ[0][0]
					else:
						happening[3]=math.inf
				break


def CpuDone():
	global clock1, subscribers_in_system, visitors_in_system, happening, cQ, dQ, oQ, bQ, subscribers_in_out, visitors_in_out

	clock1=happening[0]

	cpu["numberCompletions"] = cpu["numberCompletions"] + 1
	cpu["sumTimeLength"] = cpu["sumTimeLength"] + (clock1 - cpu["oldClock"])*len(cQ)
	cpu["sumBusyTime"] = cpu["sumBusyTime"] + (clock1 - cpu["oldClock"])*min(len(cQ),cpu["numberServers"])
	cpu["oldClock"] = clock1


	head=cQ[0]
	ReductDatTime(head[0])
	cQ.pop(0)

	if(head[5]!=0):
		if(head[4]=='S'):
			cpu["numberCompletions_s"] = cpu["numberCompletions_s"] + 1
			request = random.expovariate(1/disk["meanService_s"])
		else:
			request = random.expovariate(1/disk["meanService_v"])
		tup = [request, head[1], 1, head[3], head[4], head[5]-1]
		updateQueue(1, tup)

	else:
		cpu["numberCompletions_v"] = cpu["numberCompletions_v"] + 1
		if(head[4]=='S'):
			request = random.expovariate(1/out["meanService_s"])
			subscribers_in_system-=1
			subscribers_in_out+=1
		else:
			request = random.expovariate(1/out["meanService_v"])	
			visitors_in_system-=1
			visitors_in_out+=1	
		tup= [request, head[1], 3, head[3], head[4], head[5]]
		updateQueue(3, tup)

	if cQ:
		happening[0] = clock1 + cQ[0][0]
	else:
		happening[0]=math.inf

def DiskDone():
	global clock1, happening, cQ, dQ, oQ, bQ

	clock1=happening[1]

	disk["numberCompletions"] = disk["numberCompletions"] + 1
	disk["sumTimeLength"] = disk["sumTimeLength"] + (clock1 - disk["oldClock"])*len(dQ)
	disk["sumBusyTime"] = disk["sumBusyTime"] + (clock1 - disk["oldClock"])*min(len(dQ),disk["numberServers"])
	disk["oldClock"] = clock1

	head=dQ[0]
	ReductDatTime(head[0])
	dQ.pop(0)

	if(head[4]=='S'):
		disk["numberCompletions_s"] = disk["numberCompletions_s"] + 1
		request = -(cpu["meanService_s"]/4) * math.log(random.random()*random.random()*random.random()*random.random())
	else:
		disk["numberCompletions_v"] = disk["numberCompletions_v"] + 1
		request = -(cpu["meanService_v"]/4) * math.log(random.random()*random.random()*random.random()*random.random())

	head[0]=request
	head[2]=0
	updateQueue(0, head)
	if dQ:
		happening[1] = clock1 + dQ[0][0]
	else:
		happening[1] = math.inf

def OutDone():
	global clock1, happening, cQ, dQ, oQ, bQ, subscribers_in_out, visitors_in_out

	clock1=happening[3]

	out["numberCompletions"] = out["numberCompletions"] + 1
	out["sumTimeLength"] = out["sumTimeLength"] + (clock1 - out["oldClock"])*len(oQ)
	out["sumBusyTime"] = out["sumBusyTime"] + (clock1 - out["oldClock"])*min(len(oQ),out["numberServers"])
	out["oldClock"] = clock1

	head=oQ[0]
	ReductDatTime(head[0])

	oQ.pop(0)

	if oQ:
		happening[3] = clock1 + oQ[0][0]
	else:
		happening[3] = math.inf 
	index=head[1]

	if head[4]=='V': #get it out of balkers
		visitors_in_out-=1
		out["numberCompletions_v"] = out["numberCompletions_v"] + 1
		for i in range(len(bQ)):
			if bQ[i][1]==index:
				bQ.pop(i)
				if i==0:
					if bQ:
						happening[2]=clock1+bQ[0][0]
					else: 
						happening[2]=math.inf
				break
	else:
		subscribers_in_out-=1
		out["numberCompletions_s"] = out["numberCompletions_s"] + 1

def BalkingDone():
	global total_balkers_v, clock1, happening, cQ, dQ, oQ, bQ
	clock1=happening[2]
	total_balkers_v+=1
	head=bQ[0]

	ReductDatTime(head[0])

	pop_it(head[1])
	bQ.pop(0)

	if bQ:
		happening[2]=clock1 + bQ[0][0]
	else:
		happening[2]=math.inf


def checkCycle():
	global numberCycles, cycleLength, endCycle, timeCycleStarted, sumcl, sumclsq, clock1, happening, cQ, dQ, oQ, bQ

	endCycle=0
	if numberEvents>0 and not dQ and not cQ and not oQ:
		endCycle=1
		numberCycles+=1
		cycleLength = clock1 - timeCycleStarted;
		timeCycleStarted = clock1;
		sumcl = sumcl + cycleLength;
		sumclsq = sumclsq + (cycleLength*cycleLength);

		cpu["sumTimeLength"] = cpu["sumTimeLength"] + (clock1 - cpu["oldClock"])*len(cQ)
		cpu["sumBusyTime"] = (cpu["sumBusyTime"] + (clock1-cpu["oldClock"])*min(len(cQ),cpu["numberServers"]) )/ cpu["numberServers"] 
		cpu["oldClock"] = clock1

		cpu["bt"] = cpu["bt"] + cpu["sumBusyTime"]
		cpu["tl"] = cpu["tl"] + cpu["sumTimeLength"]
		cpu["nc"] = cpu["nc"] + cpu["numberCompletions"]
		cpu["nc_s"] = cpu["nc_s"] + cpu["numberCompletions_s"]
		cpu["nc_v"] = cpu["nc_v"] + cpu["numberCompletions_v"]
		cpu["sumBusyTime"]=0.0
		cpu["tlsq"] = cpu["tlsq"] + math.sqrt(cpu["sumTimeLength"]*cpu["sumTimeLength"])
		cpu["tlxcl"] = cpu["tlxcl"] + cpu["sumTimeLength"]*cycleLength
		cpu["tlxnc"] = cpu["tlxnc"] + cpu["sumTimeLength"]*cpu["numberCompletions"]
		cpu["ncsq"] = cpu["ncsq"] + math.sqrt(cpu["numberCompletions"])
		cpu["ncsq_s"] = cpu["ncsq"] + math.sqrt(cpu["numberCompletions_s"])
		cpu["numberCompletions"] = 0
		cpu["numberCompletions_s"] = 0
		cpu["sumTimeLength"] = 0


		disk["sumTimeLength"] = disk["sumTimeLength"] + (clock1 - disk["oldClock"])*len(dQ)
		disk["sumBusyTime"] = (disk["sumBusyTime"] + (clock1-disk["oldClock"])*min(len(dQ),disk["numberServers"]) )/ disk["numberServers"] 
		disk["oldClock"] = clock1

		disk["bt"] = disk["bt"] + disk["sumBusyTime"]
		disk["tl"] = disk["tl"] + disk["sumTimeLength"]
		disk["nc"] = disk["nc"] + disk["numberCompletions"]
		disk["nc_s"] = disk["nc_s"] + disk["numberCompletions_s"]
		disk["nc_v"] = disk["nc_v"] + disk["numberCompletions_v"]
		disk["sumBusyTime"]=0.0
		disk["tlsq"] = disk["tlsq"] + math.sqrt(disk["sumTimeLength"]*disk["sumTimeLength"])
		disk["tlxcl"] = disk["tlxcl"] + disk["sumTimeLength"]*cycleLength
		disk["tlxnc"] = disk["tlxnc"] + disk["sumTimeLength"]*disk["numberCompletions"]
		disk["ncsq"] = disk["ncsq"] + math.sqrt(disk["numberCompletions"])
		disk["ncsq_s"] = disk["ncsq"] + math.sqrt(disk["numberCompletions_s"])
		disk["numberCompletions"] = 0
		disk["sumTimeLength"] = 0

		out["sumTimeLength"] = out["sumTimeLength"] + (clock1 - out["oldClock"])*len(cQ)
		out["sumBusyTime"] = (out["sumBusyTime"] + (clock1-out["oldClock"])*min(len(cQ),out["numberServers"]) )/ out["numberServers"] 
		out["oldClock"] = clock1

		out["bt"] = out["bt"] + out["sumBusyTime"]
		out["tl"] = out["tl"] + out["sumTimeLength"]
		out["nc"] = out["nc"] + out["numberCompletions"]
		out["nc_s"] = out["nc_s"] + out["numberCompletions_s"]
		out["nc_v"] = out["nc_v"] + out["numberCompletions_v"]
		out["sumBusyTime"]=0.0
		out["tlsq"] = out["tlsq"] + math.sqrt(out["sumTimeLength"]*out["sumTimeLength"])
		out["tlxcl"] = out["tlxcl"] + out["sumTimeLength"]*cycleLength
		out["tlxnc"] = out["tlxnc"] + out["sumTimeLength"]*out["numberCompletions"]
		out["ncsq"] = out["ncsq"] + math.sqrt(out["numberCompletions"])
		out["ncsq_s"] = out["ncsq"] + math.sqrt(out["numberCompletions_s"])
		out["numberCompletions"] = 0
		out["sumTimeLength"] = 0


def addCPU():
	global total_req, subscribers_in_system, total_s, visitors_in_system, total_v, total_balkers_s, clock1, happening, cQ, dQ, oQ, bQ, N_s, N_v

	ReductDatTime(happening[4]-clock1)
	clock1=happening[4]

	cpu["sumTimeLength"] = cpu["sumTimeLength"] + (clock1 - cpu["oldClock"])*len(cQ)
	cpu["sumBusyTime"] = cpu["sumBusyTime"] + (clock1 - cpu["oldClock"])*min(len(cQ),cpu["numberServers"])
	cpu["oldClock"] = clock1

	clients = random.choices(['Subscriber', 'Visitor'], [0.7, 0.3], k=1)
	client=clients[0]

	if (client=='Subscriber'):
		request= -(cpu["meanService_s"]/4) * math.log(random.random()*random.random()*random.random()*random.random())
		type1='S'
		id1=total_req+1
		xi = random.normalvariate(40.5, 6)
		visits=  int(random.expovariate(1/meandiskvisits_s))
		#tup = [request, id, station, ktofli, type, visits_to_disk]
		tup = [request, id1, 0, xi, type1, visits]
		if(xi>=subscribers_in_system+visitors_in_system):
			updateQueue(0, tup)
			subscribers_in_system=subscribers_in_system+1
			total_s = total_s + 1 

		else:
			total_balkers_s=total_balkers_s+1;

	else: #visitor
		request = -(cpu["meanService_v"]/4) * math.log(random.random()*random.random()*random.random()*random.random())
		type1='V'
		id1=total_req+1
		xi = 1000 * 4.2 / ((1-random.random())**(1/1.5)) #also turn it to msec
		visits = int(random.expovariate(1/meandiskvisits_v) )
		tup = [request, id1, 0, xi, type1, visits]
		updateQueue(0, tup)
		updateQueue(2, tup)
		visitors_in_system=visitors_in_system+1
		total_v = total_v + 1 


	if cQ:
		happening[0]= clock1 + cQ[0][0]
	else:
		happening[0] = math.inf

	total_req+=1

	#next addition to CPU in
	happening[4] -= 1000*math.log(random.random())/2.4#np.random.exponential(2.4)*1000 #turn it to msec

N_s=0
N_v=0
subscribers_in_out=0
visitors_in_out=0
numberEvents = 0
endCycle=0
timeCycleStarted=0.0
sumcl=0.0
sumclsq=0.0
total_req=0
subscribers_in_system=0		#currently in the system
visitors_in_system=0
total_s=0	#total in system
total_v=0
total_balkers_v=0	
total_balkers_s=0	

timeCycleStarted=0.0

#first visitor
numberCycles=0
happening[4]= -1000*math.log(random.random())/2.4 #np.random.exponential(2.4)*1000 #turn it to msec
clock1=happening[4]
addCPU()


##########################################################################
########################### THE LOOP #####################################
##########################################################################

while(numberCycles<1000): 

	####################################################
	event = minf(happening)
	oldie=clock1
	if happening[event]!=math.inf:
		if(event==0):
			CpuDone()
		elif(event==1):
			DiskDone()
		elif(event==2):
			BalkingDone()
		elif (event==3):
			OutDone()
		elif (event==4):
			addCPU()
	N_s+=(subscribers_in_system+subscribers_in_out)*(clock1-oldie)
	N_v+=(visitors_in_system+visitors_in_out)*(clock1-oldie)

	numberEvents+=1
	checkCycle();
	
	if(numberCycles%20==0 and numberCycles!=0):
		cycleLength = sumcl/numberCycles;
		nocycm1 = numberCycles - 1;
		varcl = (sumclsq-math.sqrt(sumcl)/numberCycles)/nocycm1;

		#vathmos empisotshnhs 1-a=0.95, z_{1-a/2}=1.96
		if(disk["nc"]>0.0):
			#disk util
			disk["util"]=disk["bt"]/sumcl;

			disk["varnc"] = (disk["ncsq"]-math.sqrt(disk["nc_s"])/numberCycles)/nocycm1; 

			disk["ql"] = disk["tl"]/sumcl
			disk["vartl"] = (disk["tlsq"]-math.sqrt(disk["tl"])/numberCycles)/nocycm1
			disk["covartlcl"]=(disk["tlxcl"]-disk["tl"]*sumcl/ numberCycles)/nocycm1
			#disk["dql"] = 1.96*math.sqrt((disk["vartl"]-2*disk["ql"]*disk["covartlcl"] +  math.sqrt(disk["ql"])*varcl)/numberCycles)/cycleLength

			disk["qt_s"] = disk["tl"]/disk["nc_s"]
			disk["covartlnc"]=(disk["tlxnc"]-disk["tl"]*disk["nc"] / numberCycles)/nocycm1
			#disk["dqt"] = 1.96*math.sqrt((disk["vartl"]-2*disk["qt"]*disk["covartlnc"] +  math.sqrt(disk["qt"])*disk["varnc"])/numberCycles)/(disk["nc"]/numberCycles)
			

		if(cpu["nc"]>0.0):
			#cpu util
			cpu["util"]=cpu["bt"]/sumcl;

			cpu["ql"] = cpu["tl"]/sumcl
			cpu["vartl"] = (cpu["tlsq"]-math.sqrt(cpu["tl"]) /numberCycles)/nocycm1
			cpu["covartlcl"]=(cpu["tlxcl"]-cpu["tl"]*sumcl/ numberCycles)/nocycm1
			#cpu["dql"] = 1.96*math.sqrt((cpu["vartl"]-2*cpu["ql"]*cpu["covartlcl"] +  math.sqrt(cpu["ql"])*varcl)/numberCycles)/cycleLength
			
			cpu["varnc"] = (cpu["ncsq"]-math.sqrt(cpu["nc"])/numberCycles)/nocycm1; 

			cpu["qt_s"] = cpu["tl"]/cpu["nc_s"]
			cpu["covartlnc"]=(cpu["tlxnc"]-cpu["tl"]*cpu["nc"]/numberCycles)/nocycm1 
			#cpu["dqt"] = 1.96*math.sqrt((cpu["vartl"]-2*cpu["qt"]*cpu["covartlnc"] +  math.sqrt(cpu["qt"])*cpu["varnc"])/numberCycles)/(cpu["nc"]/numberCycles)
			#cpu["dqt_s"] = 1.96*math.sqrt((cpu["vartl"]-2*cpu["qt"]*cpu["covartlnc"] +  math.sqrt(cpu["qt_s"])*cpu["varnc"])/numberCycles)/(cpu["nc_s"]/numberCycles)


			if(cpu["dqt_s"] < 0.1*cpu["qt_s"]):
				break
		
		if(out["nc"]>0.0):
			#out util
			out["util"]=out["bt"]/sumcl;

			out["varnc"] = (out["ncsq"]-math.sqrt(out["nc_s"])/numberCycles)/nocycm1; 

			out["ql"] = out["tl"]/sumcl
			out["vartl"] = (out["tlsq"]-math.sqrt(out["tl"]) /numberCycles)/nocycm1
			out["covartlcl"]=(out["tlxcl"]-out["tl"]*sumcl/ numberCycles)/nocycm1
			#out["dql"] = 1.96*math.sqrt((out["vartl"]-2*out["ql"]*out["covartlcl"] +  math.sqrt(out["ql"])*varcl)/numberCycles)/cycleLength

			out["qt_s"] = out["tl"]/out["nc_s"]
			qt_s = out["tl"]/out["nc_v"]
			out["covartlnc"]=(out["tlxnc"]-out["tl"]*out["nc"]/numberCycles)/nocycm1 
			#out["dqt"] = 1.96*math.sqrt((out["vartl"]-2*out["qt"]*out["covartlnc"] +  math.sqrt(out["qt"])*out["varnc"])/numberCycles)/(out["nc"]/numberCycles)

print("Response Time for Subscriber: ", N_s/total_s, "msec")
print("Response Time for Visitor: ", N_v/total_v, "msec")
print("Subscribers that abandoned: ", 100*total_balkers_s / (total_s+total_balkers_s), "%")
print("Visitors that abandoned: ", 100*total_balkers_v / (total_v+total_balkers_v), "%")
print("CPU Utilization: ", 100*cpu["util"], "%")
print("Disk Utilization: ", 100*disk["util"], "%")
print("Outing queue Utilization: ", 100*out["util"], "%")
