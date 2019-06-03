/*
 * Title:        CloudSim Toolkit
 * Description:  CloudSim (Cloud Simulation) Toolkit for Modeling and Simulation of Clouds
 * Licence:      GPL - http://www.gnu.org/copyleft/gpl.html
 *
 * Copyright (c) 2009-2012, The University of Melbourne, Australia
 */

package org.cloudbus.cloudsim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.core.CloudSimTags;
import org.cloudbus.cloudsim.core.SimEntity;
import org.cloudbus.cloudsim.core.SimEvent;
import org.cloudbus.cloudsim.lists.CloudletList;
import org.cloudbus.cloudsim.lists.VmList;

/**
 * DatacentreBroker represents a broker acting on behalf of a user. It hides VM management, as vm
 * creation, submission of cloudlets to VMs and destruction of VMs.
 * 
 * @author Rodrigo N. Calheiros
 * @author Anton Beloglazov
 * @since CloudSim Toolkit 1.0
 */
public class DatacenterBroker extends SimEntity {

	/** The list of VMs submitted to be managed by the broker. */
	protected List<? extends Vm> vmList;

	/** The list of VMs created by the broker. */
	protected List<? extends Vm> vmsCreatedList;

	/** The list of cloudlet submitted to the broker. 
         * @see #submitCloudletList(java.util.List) 
         */
	protected List<? extends Cloudlet> cloudletList;

	/** The list of submitted cloudlets. */
	protected List<? extends Cloudlet> cloudletSubmittedList;

	/** The list of received cloudlet. */
	protected List<? extends Cloudlet> cloudletReceivedList;

	/** The number of submitted cloudlets. */
	protected int cloudletsSubmitted;

	/** The number of requests to create VM. */
	protected int vmsRequested;

	/** The number of acknowledges (ACKs) sent in response to
         * VM creation requests. */
	protected int vmsAcks;

	/** The number of destroyed VMs. */
	protected int vmsDestroyed;

	/** The id's list of available datacenters. */
	protected List<Integer> datacenterIdsList;

	/** The list of datacenters where was requested to place VMs. */
	protected List<Integer> datacenterRequestedIdsList;

	/** The vms to datacenters map, where each key is a VM id
         * and each value is the datacenter id whwere the VM is placed. */
	protected Map<Integer, Integer> vmsToDatacentersMap;

	/** The datacenter characteristics map where each key
         * is a datacenter id and each value is its characteristics.. */
	protected Map<Integer, DatacenterCharacteristics> datacenterCharacteristicsList;

	/**
	 * Created a new DatacenterBroker object.
	 * 
	 * @param name name to be associated with this entity (as required by {@link SimEntity} class)
	 * @throws Exception the exception
	 * @pre name != null
	 * @post $none
	 */
	public DatacenterBroker(String name) throws Exception {
		super(name);

		setVmList(new ArrayList<Vm>());
		setVmsCreatedList(new ArrayList<Vm>());
		setCloudletList(new ArrayList<Cloudlet>());
		setCloudletSubmittedList(new ArrayList<Cloudlet>());
		setCloudletReceivedList(new ArrayList<Cloudlet>());

		cloudletsSubmitted = 0;
		setVmsRequested(0);
		setVmsAcks(0);
		setVmsDestroyed(0);

		setDatacenterIdsList(new LinkedList<Integer>());
		setDatacenterRequestedIdsList(new ArrayList<Integer>());
		setVmsToDatacentersMap(new HashMap<Integer, Integer>());
		setDatacenterCharacteristicsList(new HashMap<Integer, DatacenterCharacteristics>());
	}

	/**
	 * This method is used to send to the broker the list with virtual machines that must be
	 * created.
	 * 
	 * @param list the list
	 * @pre list !=null
	 * @post $none
	 */
	public void submitVmList(List<? extends Vm> list) {
		getVmList().addAll(list);
	}

	/**
	 * This method is used to send to the broker the list of cloudlets.
	 * 
	 * @param list the list
	 * @pre list !=null
	 * @post $none
         * 
         * @todo The name of the method is confused with the {@link #submitCloudlets()},
         * that in fact submit cloudlets to VMs. The term "submit" is being used
         * ambiguously. The method {@link #submitCloudlets()} would be named "sendCloudletsToVMs"
         * 
         * The method {@link #submitVmList(java.util.List)} may have
         * be checked too.
	 */
	public void submitCloudletList(List<? extends Cloudlet> list) {
		getCloudletList().addAll(list);
	}

	/**
	 * Specifies that a given cloudlet must run in a specific virtual machine.
	 * 
	 * @param cloudletId ID of the cloudlet being bount to a vm
	 * @param vmId the vm id
	 * @pre cloudletId > 0
	 * @pre id > 0
	 * @post $none
	 */
	public void bindCloudletToVm(int cloudletId, int vmId) {
		CloudletList.getById(getCloudletList(), cloudletId).setVmId(vmId);
	}
	
	public void bindCloudletToVmAnt(int vmNum,int cloudletNum) {
		int iteratorNum = 50;//迭代次数
		int antNum = 100;//蚂蚁数量
		int randomAnt = 15;//随机的蚂蚁数量
		double timeMatrix[][] = new double[cloudletNum][vmNum];//处理时间矩阵
		double pheromoneMatrix[][] = new double[cloudletNum][vmNum];//信息素矩阵
		int pathMatrix[][][]= new int[antNum][cloudletNum][vmNum];//所有蚂蚁的路径
		int pathMatrix1[][][]= new int[antNum][cloudletNum][vmNum];//所有蚂蚁的路径
		double timeAllAnt[] = new double[antNum];//每只蚂蚁的处理时间
		double timeAllAnt1[] = new double[antNum];
		
		init(vmNum, cloudletNum, timeMatrix, pheromoneMatrix);//初始化
		double pheromoneMatrix1[][] = new double[cloudletNum][vmNum];
		init(vmNum, cloudletNum, timeMatrix, pheromoneMatrix1);//初始化
		Log.printLine("原时间为："+calOriginTime(timeMatrix));
		itProcess(randomAnt,iteratorNum, antNum, vmNum, cloudletNum, pathMatrix1, timeMatrix, pheromoneMatrix1, timeAllAnt1,false);//迭代过程
		Log.printLine("原时间为："+calOriginTime(timeMatrix));
		itProcess(randomAnt,iteratorNum, antNum, vmNum, cloudletNum, pathMatrix, timeMatrix, pheromoneMatrix, timeAllAnt, true);
		int best = findminIndex(timeAllAnt);
		for(int i = 0;i<cloudletNum;i++) {
			for(int j = 0;j< vmNum;j++) {
				if(pathMatrix[best][i][j]==1) {
					bindCloudletToVm(i, j);
				}
			}
		}//将找到的最短路径绑定
	}
	/**
	 * 初始化时间矩阵和信息素矩阵
	 * @param vmNum
	 * @param cloudletNum
	 * @param timeMatrix
	 * @param pheromoneMatrix
	 */
	public void init(int vmNum,int cloudletNum,double timeMatrix[][],double pheromoneMatrix[][]) {
		for(int i=0;i<cloudletNum;i++) {
			for(int j=0;j<vmNum;j++) {
				timeMatrix[i][j]=(double)CloudletList.getById(getCloudletList(), i).getCloudletLength()/vmList.get(j).getMips();
			}
		}//初始化时间矩阵
		
		for(int i=0;i<cloudletNum;i++) {
			for(int j=0;j<vmNum;j++) {
				pheromoneMatrix[i][j]=1;
			}
		}//初始化信息素矩阵
	}
	/**
	 * 计算原时间
	 * @param timeMatrix
	 * @return
	 */
	public double calOriginTime(double timeMatrix[][]) {
		double originTime = 0;//原始时间
		double nodeTime[] = new double[timeMatrix[0].length];
		int cloudletCount = 0,vmCount = 0;
		while(cloudletCount<timeMatrix.length) {
			nodeTime[vmCount] += timeMatrix[cloudletCount][vmCount];
			cloudletCount++;
			vmCount++;
			if(vmCount==timeMatrix[0].length) {
				vmCount = 0;
			}
		}
		for(double max:nodeTime) {
			if(max>originTime) {
				originTime = max;
			}
		}
		return originTime;
	}
	/**
	 * 返回信息素矩阵每行中信息素最大的下标数组
	 * @param pheromoneMatrix
	 * @return
	 */
	public int[] findMaxPheromone(double[][] pheromoneMatrix) {
		ArrayList <Integer> maxIndex = new ArrayList<>();
		double maxPheromone = -1;
		int maxPheromoneMatrix[] = new int[pheromoneMatrix.length];
		for(int i = 0;i < pheromoneMatrix.length;i++) {
			for(int j = 0;j<pheromoneMatrix[0].length;j++) {
				if(pheromoneMatrix[i][j]>maxPheromone) {
					maxPheromone = pheromoneMatrix[i][j];
					maxIndex.clear();
					maxIndex.add(j);
				}else if(pheromoneMatrix[i][j]==maxPheromone) {
					maxIndex.add(j);
				}
			}
			maxPheromoneMatrix[i] = maxIndex.get((int)(Math.random()*maxIndex.size()));
			maxIndex.clear();
			maxPheromone = -1;
		}
		return maxPheromoneMatrix;
	}
	/**
	 * 计算单次迭代中每只蚂蚁的时间
	 * @param antNum
	 * @param timeMatrix
	 * @param pathMatrix
	 * @return
	 */
	public double[] maxTime(int antNum,double timeMatrix[][],int pathMatrix[][][]){
		double timeAllAnt[] = new double[antNum];
		double maxTime = 0;
		for(int antCount = 0;antCount<antNum;antCount++) {
			for(int j = 0;j<pathMatrix[0][0].length;j++) {
				double nodeTime = 0;
				for(int i = 0;i<pathMatrix[0].length;i++) {
					if(pathMatrix[antCount][i][j]==1) {
						nodeTime += timeMatrix[i][j];
					}
				}
				if(nodeTime>maxTime)
					maxTime = nodeTime;
			}
			timeAllAnt[antCount]=maxTime;
			maxTime = 0;
		}
		return timeAllAnt;
	}
	/**
	 * 更新信息素矩阵
	 * @param pheromoneMatrix
	 * @param timeAllAnt
	 * @param pathMatrix
	 */
	public void updatePheromoneMatrix(double pheromoneMatrix[][],double timeAllAnt[],int pathMatrix[][][],double timeMatrix[][],boolean isDecline,boolean improve) {
		double loss = 0.5;
		for(int i = 0;i<pheromoneMatrix.length;i++) {
			for(int j =0;j<pheromoneMatrix[0].length;j++) {
				pheromoneMatrix[i][j] *= loss;
			}
		}
		//挥发
		for(int antCount = 0;antCount<timeAllAnt.length;antCount++) {
			ArrayList <Integer> cloudlets = new ArrayList<>();
			double[] timeArray = new double[pheromoneMatrix[0].length];
			for(int j = 0;j<pheromoneMatrix[0].length;j++) {
				for(int i = 0;i<pheromoneMatrix.length;i++) {
					if(pathMatrix[antCount][i][j]==1) {
						cloudlets.add(i);
					}
				}
				for(int cloudletsCount = 0;cloudletsCount<cloudlets.size();cloudletsCount++) {
					timeArray[j] += timeMatrix[cloudlets.get(cloudletsCount)][j];
				}
				for(int i = 0;i<pheromoneMatrix.length;i++) {
					if(pathMatrix[antCount][i][j]==1) {
						pheromoneMatrix[i][j] += 2/timeArray[j];
					}
				}
				cloudlets.clear();//虚拟机上的任务列表清空
			}
		}
		if (improve == true) {
			if (isDecline == true) {
				int minIndex = 0;
				double min = 1000;
				for (int i = 0; i < timeAllAnt.length; i++) {
					if (timeAllAnt[i] < min) {
						min = timeAllAnt[i];
						minIndex = i;
					}
				}
				for (int i = 0; i < pheromoneMatrix.length; i++) {
					for (int j = 0; j < pheromoneMatrix[0].length; j++) {
						if (pathMatrix[minIndex][i][j] == 1) {
							pheromoneMatrix[i][j] *= 1.33;
						}
					}
				}
				isDecline = false;
			}
		}
	}
	/**
	 * 找到所有蚂蚁中最短的时间
	 * @param timeAllAnt
	 * @return
	 */
	public double findMin(double timeAllAnt[]) {
		double min = 1000;
		for(double a:timeAllAnt) {
			if(a<min)min = a;
		}
		return min;
	}
	/**
	 * 找到最短路径蚂蚁的序号
	 * @param timeAllAnt
	 * @return
	 */
	public int findminIndex(double timeAllAnt[]) {
		double min = 1000;
		int minIndex = 0;
		for(int i = 0;i<timeAllAnt.length;i++) {
			if(timeAllAnt[i]<min) {
				min = timeAllAnt[i];
				minIndex = i;
			}
		}
		return minIndex;
	}
	/**
	 * 计算云任务在各个虚拟机上分配的概率
	 * @param cloudletCount
	 * @param pheromoneMatrix
	 * @param timeMatrix
	 * @return
	 */
	public double[] calChance(int cloudletCount,double pheromoneMatrix[][],double timeMatrix[][]) {
			double alpha = 3.2;
			double beta = 4.5;
			double denominator = 0;
			double chance[] = new double[pheromoneMatrix[0].length];
				for(int j = 0;j<pheromoneMatrix[0].length;j++) {
					denominator += Math.pow(pheromoneMatrix[cloudletCount][j],alpha)*Math.pow(1/timeMatrix[cloudletCount][j], beta);
				}
				for(int j = 0;j<pheromoneMatrix[0].length;j++) {
					chance[j] = Math.pow(pheromoneMatrix[cloudletCount][j],alpha)*Math.pow(1/timeMatrix[cloudletCount][j], beta)/denominator;
				}
			return chance;
	}
	/**
	 * 根据概率选择路径
	 * @param antCount
	 * @param cloudCount
	 * @param chance
	 * @param pathMatrix
	 */
	public void choosePath(int antCount, int cloudletCount, double[] chance, int pathMatrix[][][]) {
		double rand = Math.random();
		double nowchance = 0;
		int i = 0;
		while(nowchance<=rand) {
			nowchance += chance[i];
			i++;
		}
		pathMatrix[antCount][cloudletCount][i-1] = 1;
	}
	/**
	 * 迭代过程
	 * @param randomAnt
	 * @param iteratorNum
	 * @param antNum
	 * @param vmNum
	 * @param cloudletNum
	 * @param pathMatrix
	 * @param timeMatrix
	 * @param pheromoneMatrix
	 * @param timeAllAnt
	 */
	public void itProcess(int randomAnt,int iteratorNum,int antNum,int vmNum,int cloudletNum,int pathMatrix[][][],double timeMatrix[][],double pheromoneMatrix[][],double timeAllAnt[],boolean improve) {
		double lastTime = 0;//上一次迭代的最短时间
		boolean isDecline = false;//是否下降的布尔值
		
		for(int itCount=0;itCount<iteratorNum;itCount++) {
			Log.print("迭代次数："+itCount+"\n");
			for(int antCount=0;antCount<antNum-randomAnt;antCount++) {
				for(int cloudletCount = 0;cloudletCount<cloudletNum;cloudletCount++) {
					choosePath(antCount, cloudletCount, calChance(cloudletCount, pheromoneMatrix, timeMatrix), pathMatrix);
				}
			}//根据信息素启发的蚂蚁
			
			for(int antCount = antNum-randomAnt;antCount<antNum;antCount++) {
				for(int cloudletCount = 0;cloudletCount<cloudletNum;cloudletCount++) {
					pathMatrix[antCount][cloudletCount][(int)(Math.random()*vmNum)]=1;
				}
			}//随机蚂蚁
			
			timeAllAnt = maxTime(antNum, timeMatrix, pathMatrix);
			//计算时间
			Log.print("最小时间："+findMin(timeAllAnt)+"\n");
			if(lastTime>findMin(timeAllAnt)) {
				isDecline = true;
			}
			lastTime = findMin(timeAllAnt);
			updatePheromoneMatrix(pheromoneMatrix, timeAllAnt, pathMatrix,timeMatrix,isDecline,improve);
			//更新信息素
			if(itCount!=iteratorNum-1) {
			for(int i = 0;i<antNum;i++) {
				for(int j = 0;j<cloudletNum;j++) {
					for(int k = 0;k<vmNum;k++) {
						pathMatrix[i][j][k]=0;
					}
				}
			}
			}//清空路径矩阵
		}
	}
	@Override
	public void processEvent(SimEvent ev) {
		switch (ev.getTag()) {
		// Resource characteristics request
			case CloudSimTags.RESOURCE_CHARACTERISTICS_REQUEST:
				processResourceCharacteristicsRequest(ev);
				break;
			// Resource characteristics answer
			case CloudSimTags.RESOURCE_CHARACTERISTICS:
				processResourceCharacteristics(ev);
				break;
			// VM Creation answer
			case CloudSimTags.VM_CREATE_ACK:
				processVmCreate(ev);
				break;
			// A finished cloudlet returned
			case CloudSimTags.CLOUDLET_RETURN:
				processCloudletReturn(ev);
				break;
			// if the simulation finishes
			case CloudSimTags.END_OF_SIMULATION:
				shutdownEntity();
				break;
			// other unknown tags are processed by this method
			default:
				processOtherEvent(ev);
				break;
		}
	}

	/**
	 * Process the return of a request for the characteristics of a Datacenter.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != $null
	 * @post $none
	 */
	protected void processResourceCharacteristics(SimEvent ev) {
		DatacenterCharacteristics characteristics = (DatacenterCharacteristics) ev.getData();
		getDatacenterCharacteristicsList().put(characteristics.getId(), characteristics);

		if (getDatacenterCharacteristicsList().size() == getDatacenterIdsList().size()) {
			setDatacenterRequestedIdsList(new ArrayList<Integer>());
			createVmsInDatacenter(getDatacenterIdsList().get(0));
		}
	}

	/**
	 * Process a request for the characteristics of a PowerDatacenter.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != $null
	 * @post $none
	 */
	protected void processResourceCharacteristicsRequest(SimEvent ev) {
		setDatacenterIdsList(CloudSim.getCloudResourceList());
		setDatacenterCharacteristicsList(new HashMap<Integer, DatacenterCharacteristics>());

//		Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Cloud Resource List received with ",
//				getDatacenterIdsList().size(), " resource(s)");

		for (Integer datacenterId : getDatacenterIdsList()) {
			sendNow(datacenterId, CloudSimTags.RESOURCE_CHARACTERISTICS, getId());
		}
	}

	/**
	 * Process the ack received due to a request for VM creation.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != null
	 * @post $none
	 */
	protected void processVmCreate(SimEvent ev) {
		int[] data = (int[]) ev.getData();
		int datacenterId = data[0];
		int vmId = data[1];
		int result = data[2];

		if (result == CloudSimTags.TRUE) {
			getVmsToDatacentersMap().put(vmId, datacenterId);
			getVmsCreatedList().add(VmList.getById(getVmList(), vmId));
//			Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": VM #", vmId,
//					" has been created in Datacenter #", datacenterId, ", Host #",
//					VmList.getById(getVmsCreatedList(), vmId).getHost().getId());
		} else {
//			Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Creation of VM #", vmId,
//					" failed in Datacenter #", datacenterId);
		}

		incrementVmsAcks();

		// all the requested VMs have been created
		if (getVmsCreatedList().size() == getVmList().size() - getVmsDestroyed()) {
			submitCloudlets();
		} else {
			// all the acks received, but some VMs were not created
			if (getVmsRequested() == getVmsAcks()) {
				// find id of the next datacenter that has not been tried
				for (int nextDatacenterId : getDatacenterIdsList()) {
					if (!getDatacenterRequestedIdsList().contains(nextDatacenterId)) {
						createVmsInDatacenter(nextDatacenterId);
						return;
					}
				}

				// all datacenters already queried
				if (getVmsCreatedList().size() > 0) { // if some vm were created
					submitCloudlets();
				} else { // no vms created. abort
					Log.printLine(CloudSim.clock() + ": " + getName()
							+ ": none of the required VMs could be created. Aborting");
					finishExecution();
				}
			}
		}
	}

	/**
	 * Process a cloudlet return event.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != $null
	 * @post $none
	 */
	protected void processCloudletReturn(SimEvent ev) {
		Cloudlet cloudlet = (Cloudlet) ev.getData();
		getCloudletReceivedList().add(cloudlet);
//		Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Cloudlet ", cloudlet.getCloudletId(),
//				" received");
		cloudletsSubmitted--;
		if (getCloudletList().size() == 0 && cloudletsSubmitted == 0) { // all cloudlets executed
			Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": All Cloudlets executed. Finishing...");
			clearDatacenters();
			finishExecution();
		} else { // some cloudlets haven't finished yet
			if (getCloudletList().size() > 0 && cloudletsSubmitted == 0) {
				// all the cloudlets sent finished. It means that some bount
				// cloudlet is waiting its VM be created
				clearDatacenters();
				createVmsInDatacenter(0);
			}

		}
	}

	/**
	 * Process non-default received events that aren't processed by
         * the {@link #processEvent(org.cloudbus.cloudsim.core.SimEvent)} method.
         * This method should be overridden by subclasses in other to process
         * new defined events.
	 * 
	 * @param ev a SimEvent object
	 * @pre ev != null
	 * @post $none
         * @todo to ensure the method will be overridden, it should be defined 
         * as abstract in a super class from where new brokers have to be extended.
	 */
	protected void processOtherEvent(SimEvent ev) {
		if (ev == null) {
			Log.printConcatLine(getName(), ".processOtherEvent(): ", "Error - an event is null.");
			return;
		}

		Log.printConcatLine(getName(), ".processOtherEvent(): Error - event unknown by this DatacenterBroker.");
	}

	/**
	 * Create the submitted virtual machines in a datacenter.
	 * 
	 * @param datacenterId Id of the chosen Datacenter
	 * @pre $none
	 * @post $none
         * @see #submitVmList(java.util.List) 
	 */
	protected void createVmsInDatacenter(int datacenterId) {
		// send as much vms as possible for this datacenter before trying the next one
		int requestedVms = 0;
//		String datacenterName = CloudSim.getEntityName(datacenterId);
		for (Vm vm : getVmList()) {
			if (!getVmsToDatacentersMap().containsKey(vm.getId())) {
//				Log.printLine(CloudSim.clock() + ": " + getName() + ": Trying to Create VM #" + vm.getId()
//						+ " in " + datacenterName);
				sendNow(datacenterId, CloudSimTags.VM_CREATE_ACK, vm);
				requestedVms++;
			}
		}

		getDatacenterRequestedIdsList().add(datacenterId);

		setVmsRequested(requestedVms);
		setVmsAcks(0);
	}

	/**
	 * Submit cloudlets to the created VMs.
	 * 
	 * @pre $none
	 * @post $none
         * @see #submitCloudletList(java.util.List) 
	 */
	protected void submitCloudlets() {
		int vmIndex = 0;
		List<Cloudlet> successfullySubmitted = new ArrayList<Cloudlet>();
		for (Cloudlet cloudlet : getCloudletList()) {
			Vm vm;
			// if user didn't bind this cloudlet and it has not been executed yet
			if (cloudlet.getVmId() == -1) {
				vm = getVmsCreatedList().get(vmIndex);
			} else { // submit to the specific vm
				vm = VmList.getById(getVmsCreatedList(), cloudlet.getVmId());
				if (vm == null) { // vm was not created
					if(!Log.isDisabled()) {				    
					    Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Postponing execution of cloudlet ",
							cloudlet.getCloudletId(), ": bount VM not available");
					}
					continue;
				}
			}

			if (!Log.isDisabled()) {
//			    Log.printConcatLine(CloudSim.clock(), ": ", getName(), ": Sending cloudlet ",
//					cloudlet.getCloudletId(), " to VM #", vm.getId());
			}
			
			cloudlet.setVmId(vm.getId());
			sendNow(getVmsToDatacentersMap().get(vm.getId()), CloudSimTags.CLOUDLET_SUBMIT, cloudlet);
			cloudletsSubmitted++;
			vmIndex = (vmIndex + 1) % getVmsCreatedList().size();
			getCloudletSubmittedList().add(cloudlet);
			successfullySubmitted.add(cloudlet);
		}

		// remove submitted cloudlets from waiting list
		getCloudletList().removeAll(successfullySubmitted);
	}

	/**
	 * Destroy all virtual machines running in datacenters.
	 * 
	 * @pre $none
	 * @post $none
	 */
	protected void clearDatacenters() {
//		for (Vm vm : getVmsCreatedList()) {
//			Log.printConcatLine(CloudSim.clock(), ": " + getName(), ": Destroying VM #", vm.getId());
//			sendNow(getVmsToDatacentersMap().get(vm.getId()), CloudSimTags.VM_DESTROY, vm);
//		}

		getVmsCreatedList().clear();
	}

	/**
	 * Send an internal event communicating the end of the simulation.
	 * 
	 * @pre $none
	 * @post $none
	 */
	protected void finishExecution() {
		sendNow(getId(), CloudSimTags.END_OF_SIMULATION);
	}

	@Override
	public void shutdownEntity() {
		Log.printConcatLine(getName(), " is shutting down...");
	}

	@Override
	public void startEntity() {
		Log.printConcatLine(getName(), " is starting...");
		schedule(getId(), 0, CloudSimTags.RESOURCE_CHARACTERISTICS_REQUEST);
	}

	/**
	 * Gets the vm list.
	 * 
	 * @param <T> the generic type
	 * @return the vm list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Vm> List<T> getVmList() {
		return (List<T>) vmList;
	}

	/**
	 * Sets the vm list.
	 * 
	 * @param <T> the generic type
	 * @param vmList the new vm list
	 */
	protected <T extends Vm> void setVmList(List<T> vmList) {
		this.vmList = vmList;
	}

	/**
	 * Gets the cloudlet list.
	 * 
	 * @param <T> the generic type
	 * @return the cloudlet list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Cloudlet> List<T> getCloudletList() {
		return (List<T>) cloudletList;
	}

	/**
	 * Sets the cloudlet list.
	 * 
	 * @param <T> the generic type
	 * @param cloudletList the new cloudlet list
	 */
	protected <T extends Cloudlet> void setCloudletList(List<T> cloudletList) {
		this.cloudletList = cloudletList;
	}

	/**
	 * Gets the cloudlet submitted list.
	 * 
	 * @param <T> the generic type
	 * @return the cloudlet submitted list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Cloudlet> List<T> getCloudletSubmittedList() {
		return (List<T>) cloudletSubmittedList;
	}

	/**
	 * Sets the cloudlet submitted list.
	 * 
	 * @param <T> the generic type
	 * @param cloudletSubmittedList the new cloudlet submitted list
	 */
	protected <T extends Cloudlet> void setCloudletSubmittedList(List<T> cloudletSubmittedList) {
		this.cloudletSubmittedList = cloudletSubmittedList;
	}

	/**
	 * Gets the cloudlet received list.
	 * 
	 * @param <T> the generic type
	 * @return the cloudlet received list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Cloudlet> List<T> getCloudletReceivedList() {
		return (List<T>) cloudletReceivedList;
	}

	/**
	 * Sets the cloudlet received list.
	 * 
	 * @param <T> the generic type
	 * @param cloudletReceivedList the new cloudlet received list
	 */
	protected <T extends Cloudlet> void setCloudletReceivedList(List<T> cloudletReceivedList) {
		this.cloudletReceivedList = cloudletReceivedList;
	}

	/**
	 * Gets the vm list.
	 * 
	 * @param <T> the generic type
	 * @return the vm list
	 */
	@SuppressWarnings("unchecked")
	public <T extends Vm> List<T> getVmsCreatedList() {
		return (List<T>) vmsCreatedList;
	}

	/**
	 * Sets the vm list.
	 * 
	 * @param <T> the generic type
	 * @param vmsCreatedList the vms created list
	 */
	protected <T extends Vm> void setVmsCreatedList(List<T> vmsCreatedList) {
		this.vmsCreatedList = vmsCreatedList;
	}

	/**
	 * Gets the vms requested.
	 * 
	 * @return the vms requested
	 */
	protected int getVmsRequested() {
		return vmsRequested;
	}

	/**
	 * Sets the vms requested.
	 * 
	 * @param vmsRequested the new vms requested
	 */
	protected void setVmsRequested(int vmsRequested) {
		this.vmsRequested = vmsRequested;
	}

	/**
	 * Gets the vms acks.
	 * 
	 * @return the vms acks
	 */
	protected int getVmsAcks() {
		return vmsAcks;
	}

	/**
	 * Sets the vms acks.
	 * 
	 * @param vmsAcks the new vms acks
	 */
	protected void setVmsAcks(int vmsAcks) {
		this.vmsAcks = vmsAcks;
	}

	/**
	 * Increment the number of acknowledges (ACKs) sent in response
         * to requests of VM creation.
	 */
	protected void incrementVmsAcks() {
		vmsAcks++;
	}

	/**
	 * Gets the vms destroyed.
	 * 
	 * @return the vms destroyed
	 */
	protected int getVmsDestroyed() {
		return vmsDestroyed;
	}

	/**
	 * Sets the vms destroyed.
	 * 
	 * @param vmsDestroyed the new vms destroyed
	 */
	protected void setVmsDestroyed(int vmsDestroyed) {
		this.vmsDestroyed = vmsDestroyed;
	}

	/**
	 * Gets the datacenter ids list.
	 * 
	 * @return the datacenter ids list
	 */
	protected List<Integer> getDatacenterIdsList() {
		return datacenterIdsList;
	}

	/**
	 * Sets the datacenter ids list.
	 * 
	 * @param datacenterIdsList the new datacenter ids list
	 */
	protected void setDatacenterIdsList(List<Integer> datacenterIdsList) {
		this.datacenterIdsList = datacenterIdsList;
	}

	/**
	 * Gets the vms to datacenters map.
	 * 
	 * @return the vms to datacenters map
	 */
	protected Map<Integer, Integer> getVmsToDatacentersMap() {
		return vmsToDatacentersMap;
	}

	/**
	 * Sets the vms to datacenters map.
	 * 
	 * @param vmsToDatacentersMap the vms to datacenters map
	 */
	protected void setVmsToDatacentersMap(Map<Integer, Integer> vmsToDatacentersMap) {
		this.vmsToDatacentersMap = vmsToDatacentersMap;
	}

	/**
	 * Gets the datacenter characteristics list.
	 * 
	 * @return the datacenter characteristics list
	 */
	protected Map<Integer, DatacenterCharacteristics> getDatacenterCharacteristicsList() {
		return datacenterCharacteristicsList;
	}

	/**
	 * Sets the datacenter characteristics list.
	 * 
	 * @param datacenterCharacteristicsList the datacenter characteristics list
	 */
	protected void setDatacenterCharacteristicsList(
			Map<Integer, DatacenterCharacteristics> datacenterCharacteristicsList) {
		this.datacenterCharacteristicsList = datacenterCharacteristicsList;
	}

	/**
	 * Gets the datacenter requested ids list.
	 * 
	 * @return the datacenter requested ids list
	 */
	protected List<Integer> getDatacenterRequestedIdsList() {
		return datacenterRequestedIdsList;
	}

	/**
	 * Sets the datacenter requested ids list.
	 * 
	 * @param datacenterRequestedIdsList the new datacenter requested ids list
	 */
	protected void setDatacenterRequestedIdsList(List<Integer> datacenterRequestedIdsList) {
		this.datacenterRequestedIdsList = datacenterRequestedIdsList;
	}

}
