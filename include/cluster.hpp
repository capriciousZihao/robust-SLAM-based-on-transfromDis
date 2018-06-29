// RRR - Robust Loop Closing over Time
// Copyright (C) 2014 Y.Latif, C.Cadena, J.Neira


#ifndef CLUSTER_HPP_
#define CLUSTER_HPP_

#include <iostream>
#include <sstream>
#include <fstream> 
#include <string.h>
#include <vector>
#include <list>
#include <cstdlib>
#include <Eigen/Dense>
#include <math.h>
#include <unistd.h>
#include <algorithm>  

#include "g2o/core/eigen_types.h"
#include "g2o/types/slam2d/se2.h"
#include "utils.hpp"

using namespace Eigen; 
using namespace std;
using namespace utils;

struct cluster
{
	int startLow, startHigh;
	int endLow, endHigh, nearestLCID=-1;
	int size, id_cluster;
	std::vector<std::array<double,6>> positionserial;
	std::vector<std::vector<std::pair<int, double> > > consistentGroup;
	std::vector<double> dis_cluster_backup;
	std::vector<double> dis_new_LC;
	std::vector<double> dis_cluster_start_end;
	std::vector<std::pair<int,int > > uncert_pair_group;

	cluster(): startLow(-1), startHigh(-1), endLow(-1), endHigh(-1), size(0) {}
	cluster(int start, int end, int idofcluster) : startLow(start), startHigh(start), endLow(end), endHigh(end), size(1), id_cluster(idofcluster){}

	bool contains(int start, int end, int threshold)
	{
		return
				(
					std::abs(start-startHigh) < threshold or	std::abs(start-startLow)  < threshold
				)
				and
				(
					std::abs(end-endHigh) < threshold or std::abs(end-endLow) < threshold
				);

	}

	//return  the ID of nearest LC and the distance
	std::array<double,2> cal_distance( std::array<double,6> LoopPosition)
	{
		
		double startDIS,endDIS;
		int IDofNearestLC=0, IDDynamic=0;
		std::array<double,6> lastLI;
		std::array<double,2> ret;

		dis_new_LC.clear();

		if(positionserial.size()==0)
			cout<<"can't calculate distance because positionserial has no information"<<endl;

		for(std::vector<std::array<double,6>>::const_iterator it = positionserial.begin(),
			lendVertex = positionserial.end();it!=lendVertex;it++)
		{
			lastLI = *it;
			startDIS = sqrt(pow(lastLI[1]-LoopPosition[1], 2)+pow(lastLI[2]-LoopPosition[2], 2));
			endDIS = sqrt(pow(lastLI[4]-LoopPosition[4], 2)+pow(lastLI[5]-LoopPosition[5], 2));
			if(it == positionserial.begin())
			{
				// ret[0] = startDIS;
				ret[0] = IDofNearestLC;
 				ret[1] = endDIS+startDIS;
			}
			else if (endDIS+startDIS < ret[1])
			{
				// ret[0] = startDIS;
				ret[0] = IDofNearestLC;
				ret[1] = endDIS+startDIS;
			}
			IDofNearestLC++;
			//add distance to list
			dis_new_LC.push_back(endDIS+startDIS);
		}
		// cout<<"size of dis_new_LC: "<<dis_new_LC.size()<<endl;
		nearestLCID = ret[0];
		return ret;
	}


	int  update_distance(int idnearest, double distance)
	{
		std::array<double,6> lastLI;
		double startDIS,endDIS;
		std::array<double,2> ret;
		//there should be element in positionserial, if no, it is abnormal,so exit
		if(positionserial.size() < 1)
		{
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		//check if idnearest is within the size of dis_new_LC 
		if(dis_new_LC.size() < idnearest)
		{
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		//check if the numbers of elements in posserial\dis_cluster\dis_new_lc are consistent 
		if(dis_new_LC.size() != positionserial.size())
		{
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(dis_cluster_start_end.size() != positionserial.size()-1)
		{
			cout<<"dis_cluster_start_end size does not equal to positionserial"<<endl;
			cout<<"dis_cluster_start_end.size:"<<dis_cluster_start_end.size()<<endl;
			cout<<"positionserial.size:"<<positionserial.size()<<endl;	
			cout<<"id of current cluster: "<<id_cluster<<endl;	
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);	
			exit(0);
		}
		//make sure the direct and indirect reference of nearest distance is identical
		if(dis_new_LC[idnearest] != distance)
		{
			cout<<"new distance serial disorder!"<<endl;
			cout<<"id of current cluster: "<<id_cluster<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		//clear the backup dis vector
		dis_cluster_backup.clear();
		//if no element in dis_cluster_start_end, save the distance to it and exit
		if(dis_cluster_start_end.size() == 0)
		{
			if(positionserial.size() != 1)
			{
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			dis_cluster_backup.push_back(distance);
			return 1;
		}
		//copy dis_cluster_start_end to dis_cluster_backup
		dis_cluster_backup.assign(dis_cluster_start_end.begin(), dis_cluster_start_end.end());  

		//insert the distance into dis_cluster_backup
		if(idnearest == 0)//if the nearest lc is the start lc
		{
			// if(dis_new_LC.size() == 1)
			// 	return 
			if(dis_cluster_backup[0]<=dis_new_LC[idnearest+1])//instert to the first place
			{
				dis_cluster_backup.insert(dis_cluster_backup.begin(),dis_new_LC[idnearest]);
				return 0;
			}
			else
			{
				dis_cluster_backup.insert(dis_cluster_backup.begin(),dis_new_LC[idnearest]);
				dis_cluster_backup[1] = dis_new_LC[idnearest+1];
				return 1;
			}
		}
		else if(idnearest == (positionserial.size() - 1))//if the nearest lc is the last lc
		{
			if(idnearest != dis_cluster_backup.size())
			{
				cout<<"positionserial.size: "<<positionserial.size()<<endl;
				cout<<"dis_cluster_backup.size: "<<dis_cluster_backup.size()<<endl;
				
				cout<<"dis_cluster_start_end.size: "<<dis_cluster_start_end.size()<<endl;
				cout<<"idnearest: "<<idnearest<<endl;
				cout<<"last ID does not match the last backup distance!"<<endl;
				cout<<"id of current cluster: "<<id_cluster<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			if(dis_cluster_backup[idnearest-1]<=dis_new_LC[ idnearest-1])//instert to the last place
			{
				dis_cluster_backup.push_back(distance);
				return -1;
			}
			else
			{
				*(dis_cluster_backup.end()-1) = dis_new_LC[idnearest-1];
				dis_cluster_backup.push_back(dis_new_LC[idnearest]);
				return -2;
			}
		}
		else
		{
			if(dis_cluster_backup[idnearest-1]<=dis_new_LC[ idnearest-1])//insert after:dis be i,posserial i+1
			{
				//need to find out the effect of insert(position,val) 
				dis_cluster_backup.insert(dis_cluster_backup.begin()+idnearest, dis_new_LC[ idnearest]);
				*(dis_cluster_backup.begin()+idnearest+1) = dis_new_LC[ idnearest+1];
				return idnearest+1;
			}
			else//insert before: dis be i,posserial i
			{
				dis_cluster_backup.insert(dis_cluster_backup.begin()+idnearest, dis_new_LC[ idnearest]);
				*(dis_cluster_backup.begin()+idnearest-1) = dis_new_LC[ idnearest-1];
				return idnearest;
			}
		}

		if(dis_cluster_start_end.size()!=positionserial.size()-1)
		{
			cout<<"update distance abnormal: dis size does not match positionserial size"<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}	
	}

	void updateAfterChi2OK(std::array<double,6> Loopinfo, int insertPOsition)
	{
		double endDIS, startDIS;
		if(insertPOsition >= 0)
			positionserial.insert(positionserial.begin()+insertPOsition,Loopinfo);//update posserial
		else
		{
			if(insertPOsition == -1)
				positionserial.push_back(Loopinfo);
			else if(insertPOsition == -2)
				positionserial.insert(positionserial.end()-1,Loopinfo);
			else
			{
				  printf("This fake error is in %s on line %d\n",         __FILE__, __LINE__);
				  exit(0);
			}
				
		}

		//clear the backup dis vector
		if(positionserial.size() <= 1)
		{
			cout<<"should not add distance to :"<<id_cluster<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		dis_cluster_start_end.clear();

		dis_cluster_start_end.assign(dis_cluster_backup.begin(), dis_cluster_backup.end());  

		int x=0;
		for(std::vector<std::array<double,6>>::const_iterator it = positionserial.begin(),
			lendVertex = positionserial.end()-1;it!=lendVertex;it++)
		{
			// cout<<"debug updateAfterChi2OK x: "<<x<<endl;
			std::array<double,6> lastLI = *it, LoopPosition = *(it+1);
			startDIS = sqrt(pow(lastLI[1]-LoopPosition[1], 2)+pow(lastLI[2]-LoopPosition[2],2));
			endDIS = sqrt(pow(lastLI[4]-LoopPosition[4],2)+pow(lastLI[5]-LoopPosition[5],2));
			if (abs(endDIS+startDIS - dis_cluster_start_end[x]) > 0.001 )
			{
				cout<<"update distance abnormal: dis size does not match positionserial size"<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			x++;
		}

	}
};

class Clusterizer
{
	typedef IntPairIDMap 			LoopToClusterIDMap;//loop is a pair of Vertex ID, this map is a int pair,LC, to a int, cluster ID 
	typedef IDintPairSetMap 		ClusterIDtoLoopsMap;
	

	std::vector<std::vector< std::vector<int> > >  node_sequence;
	std::vector<std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > > transSequence_whole;
	std::vector<std::array<double,4>> VertexInf;
	std::vector<std::array<double,8>> OdoInf;
	std::map<std::pair<int, int>, std::array<double,6>> LC_Inf;	
	std::vector<double> distance;
	std::vector<cluster> _clustersFound;
	ClusterIDtoLoopsMap	clusterIDtoLoopsMap;
	LoopToClusterIDMap  loopToClusterIDMap;
	std::vector<std::vector<int> > merge_vector;
	std::vector<int> killed_cluster;
	Matrix3d displayCov;

public:
	string     nameofclusterfile;


	// Assumes that in the vector the values are passed as (start_1,end_1), (start_2,end_2), ...
	int getClusterID(const IntPair& loop)//input a pair, lC, output the ID of the cluster to which the pair belongs
	{
		return loopToClusterIDMap[loop];
	}

	bool get_cluster_node_scequence(std::vector<std::vector< std::vector<int> > > & node_sequence, const std::vector<cluster> & clustersFound)
	{

		std::vector<int> nodes_dout, nodes_din;
		std::vector<std::vector<int> > element;
		bool retu = 0;
		node_sequence.clear();
		for(int i = 0; i < clustersFound.size(); ++i)
		{
			nodes_dout.clear();
			nodes_din.clear();	
			element.clear();	
			for(int j = 0; j < clustersFound[i].positionserial.size(); ++j)
			{
				nodes_dout.push_back(clustersFound[i].positionserial[j][0]);
				nodes_din.push_back(clustersFound[i].positionserial[j][3]);
			}
			element.push_back(nodes_dout);
			element.push_back(nodes_din);	
			node_sequence.push_back(element);
			// cout<<"i: "<<i<<endl;
			
			// cout<<"clustersFound[i].positionserial.size: "<<clustersFound[i].positionserial.size()<<endl;
			
			// cout<<"nodes_dout.size: "<<nodes_dout.size()<<endl;
			// cout<<"node_sequence[i][0].size: "<<node_sequence[i][0].size()<<endl;
			// cout<<"node_sequence[i][1].size: "<<node_sequence[i][1].size()<<endl;


			// std::cout << "nodes_dout sequence of cluster["<<i<<"]:"<<endl;
			// for (std::vector< int>::iterator it=nodes_dout.begin(); it!=nodes_dout.end(); ++it)
			// 	cout <<*it<<" ";
			// cout<<endl;

			// std::cout << "nodes_din  sequence of cluster["<<i<<"]:"<<endl;
			// for (std::vector< int>::iterator it=nodes_din.begin(); it!=nodes_din.end(); ++it)
			// 	cout <<*it<<" ";
			// cout<<endl;
		}
		// exit(0);
		retu = 1;
		return retu;
	}
	bool get_cluster_node_scequence(std::vector<std::vector< std::vector<int> > > & node_sequence, const std::vector<cluster> & clustersFound,
		std::vector<std::vector<int> > & merge_vector)
	{

		std::vector<int> nodes_dout, nodes_din;
		std::vector<std::vector<int> > element;
		bool retu = 0;
		node_sequence.clear();
		for(int i = 0; i < clustersFound.size(); ++i)
		{
			nodes_dout.clear();
			nodes_din.clear();	
			element.clear();	
			for(int k =0; k<clustersFound[i].positionserial.size(); k++)
			{
				nodes_dout.push_back(clustersFound[i].positionserial[k][0]);
				nodes_din.push_back(clustersFound[i].positionserial[k][3]);
			}

			element.push_back(nodes_dout);
			element.push_back(nodes_din);	
			node_sequence.push_back(element);
			// cout<<"i: "<<i<<endl;
			
			// cout<<"clustersFound[i].positionserial.size: "<<clustersFound[i].positionserial.size()<<endl;
			
			// cout<<"nodes_dout.size: "<<nodes_dout.size()<<endl;
			// cout<<"node_sequence[i][0].size: "<<node_sequence[i][0].size()<<endl;
			// cout<<"node_sequence[i][1].size: "<<node_sequence[i][1].size()<<endl;


			// std::cout << "nodes_dout sequence of cluster["<<i<<"]:"<<endl;
			// for (std::vector< int>::iterator it=nodes_dout.begin(); it!=nodes_dout.end(); ++it)
			// 	cout <<*it<<" ";
			// cout<<endl;

			// std::cout << "nodes_din  sequence of cluster["<<i<<"]:"<<endl;
			// for (std::vector< int>::iterator it=nodes_din.begin(); it!=nodes_din.end(); ++it)
			// 	cout <<*it<<" ";
			// cout<<endl;
		}
		// exit(0);
		retu = 1;
		return retu;
	}
	// check merge 
	void detect_cluster_merge( std::vector<cluster> & clusterset, std::vector<std::vector<int > > & consistent_pair_clusterr, 
		std::vector<std::vector<std::pair<double, int> > > & nearest_dis_sequence, double thres)
	{
		std::vector<std::array<double,6>> ends_lc_cluster;
		std::vector<int> ele_consistent_pair_clusterr;
		std::array<double,2> return_cal_dis, return_chi2_test;
		std::pair<double, int> dis_ele;
		std::vector<std::pair<double, int> > sequence_dis_interClusters;
		
		consistent_pair_clusterr.clear();

		int middle_lc_id;
		//if the number of clusters <2, there is no need to merge
		if(clusterset.size()<2)
		{
			return ;
		}
		nearest_dis_sequence.clear();
		for(int i = 0;i<clusterset.size();i++)
		{
			int j = 0;
			sequence_dis_interClusters.clear();
			
			for(j = 0;j<clusterset.size();j++)
			{
				if(i == j)
				{
					dis_ele.first = 0;
					dis_ele.second = i;
				}
				else
				{
					if(clusterset[i].positionserial.size() % 2 == 0)
						middle_lc_id = (clusterset[i].positionserial.size() / 2) - 1;
					else
						middle_lc_id = ((clusterset[i].positionserial.size() + 1) / 2)-1;
					return_cal_dis = clusterset[j].cal_distance(clusterset[i].positionserial[middle_lc_id]);
					dis_ele.first = return_cal_dis[1];
					dis_ele.second = j;
				}
				sequence_dis_interClusters.push_back(dis_ele);
			}
			// // print out content:
			// std::cout << "before sort:";
			// for (std::vector<std::pair<double, int> >::iterator it=sequence_dis_interClusters.begin(); it!=sequence_dis_interClusters.end(); ++it)
			// 	cout <<"dis: "<< (*it).first<<" ID: "<<(*it).second<<"\n"<<endl;
			// // std::cout << '\n';
			std::sort(sequence_dis_interClusters.begin(), sequence_dis_interClusters.end(), cmp);
			// // print out content:
			// std::cout << "after sort:";
			// for (std::vector<std::pair<double, int> >::iterator it=sequence_dis_interClusters.begin(); it!=sequence_dis_interClusters.end(); ++it)
			// 	cout <<"dis: "<< (*it).first<<" ID: "<<(*it).second<<"\n"<<endl;
			// cout<<"i: "<<i<<" j: "<<j<<endl;
			// exit(0);
			nearest_dis_sequence.push_back(sequence_dis_interClusters);
		}
		for(int i = 0; i < clusterset.size(); i++)
		{
			if(nearest_dis_sequence[i].size() != clusterset.size())
			{
				cout<<"clusterset.size(): "<<clusterset.size()<<" nearest_dis_sequence[i].size(): "<<nearest_dis_sequence[i].size()<<" i: "<<i<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}

		for(int i = 0;i<clusterset.size();i++)
		{
			ele_consistent_pair_clusterr.clear();
			ends_lc_cluster.clear();
			ends_lc_cluster.push_back(clusterset[i].positionserial[0]);
			ends_lc_cluster.push_back(clusterset[i].positionserial.back());
			bool bit =0;
			
			for(int j = 1;j<nearest_dis_sequence[i].size();j++)
			{
				// if(j == i)
				// 	continue;		
				// ends_lc_cluster.push_back(clusterset[j].positionserial[0]);
				// ends_lc_cluster.push_back(clusterset[j].positionserial.back());
				if(ends_lc_cluster.size() != 2)
				{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
				if(clusterset[nearest_dis_sequence[i][j].second].positionserial.size() == 1)
					continue;
				for(int k = 0; k < 2; k++)
				{
					if(k<2)
					{
						//calculate distance ,update distacne ,chi2 test

						return_cal_dis = clusterset[nearest_dis_sequence[i][j].second].cal_distance(ends_lc_cluster[k]);

						clusterset[nearest_dis_sequence[i][j].second].update_distance(round(return_cal_dis[0]), return_cal_dis[1]);

						if(clusterset[nearest_dis_sequence[i][j].second].dis_cluster_backup.size()-clusterset[nearest_dis_sequence[i][j].second].dis_cluster_start_end.size() != 1)
						{
							cout<<"i: "<<i<<" j: "<<j<<" dis backup size: "<<clusterset[j].dis_cluster_backup.size()<<
								"  dis cluster start end: "<<clusterset[j].dis_cluster_start_end.size()<<endl;
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
						}
						return_chi2_test = chi2_test(clusterset[nearest_dis_sequence[i][j].second].dis_cluster_backup, 
							clusterset[nearest_dis_sequence[i][j].second].dis_cluster_start_end);

						double sum = std::accumulate(std::begin(clusterset[nearest_dis_sequence[i][j].second].dis_cluster_backup), 
							std::end(clusterset[nearest_dis_sequence[i][j].second].dis_cluster_backup), 0.0);  
						double mean =  sum / clusterset[nearest_dis_sequence[i][j].second].dis_cluster_backup.size(); //均值  
						double single_chi_statis = (return_cal_dis[1]-mean)*(return_cal_dis[1]-mean)/mean;
						// if((return_chi2_test[0] > 0.05 and abs(return_chi2_test[0] - return_chi2_test[1]) < thres) or (return_cal_dis[1] < 2))
						if((single_chi_statis < 3.84) or (return_cal_dis[1] < 2))
						{
							if((clusterset[i].positionserial.size() == 2) and (k == 0))
							{
								bit =1;
								continue;
							}
							if((clusterset[i].positionserial.size() == 2) and (bit == 1))
							{
								ele_consistent_pair_clusterr.push_back(nearest_dis_sequence[i][j].second);
								bit = 0;
								break;
							}
							ele_consistent_pair_clusterr.push_back(nearest_dis_sequence[i][j].second);
							break;
						}
					}
				}

			}
			consistent_pair_clusterr.push_back(ele_consistent_pair_clusterr);
		}
		for(int i = 0; i < consistent_pair_clusterr.size(); i++)
		{
			// if(consistent_pair_clusterr[i].size() == 0)
			// {
			// 	if()
			// 	{

			// 	}
			// 	for(int j = 0; j < consistent_pair_clusterr.size(); j++)
			// 	{
			// 		if(consistent_pair_clusterr[j].size() == 0)
			// 			continue;
			// 		//find i  
			// 		vector<int>::iterator iter=find(consistent_pair_clusterr[j].begin(),consistent_pair_clusterr[j].end(),i);  
					      
			// 		//delete i 
			// 		if(iter!=consistent_pair_clusterr[j].end())
			// 			consistent_pair_clusterr[j].erase(iter);  
			// 	}
			// }

		}


		// printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
	}

	void setClusterID(const IntPair& loopClosure, int ID)
	{
		if(loopToClusterIDMap.find(loopClosure)!=loopToClusterIDMap.end())
		{
			int oldID = loopToClusterIDMap[loopClosure];

//			for(IntPairSet::iterator
//					it = clusterIDtoLoopsMap[oldID].begin(),
//					end = clusterIDtoLoopsMap[oldID].end();
//					it!=end;
//					it++)
//
//			{
				//if (*it == loopClosure){
					clusterIDtoLoopsMap[oldID].erase(loopClosure);
				//	loopToClusterIDMap[*it] = ID;
				//	break;
				//}
//			}
			clusterIDtoLoopsMap[ID].insert(loopClosure);
			loopToClusterIDMap[loopClosure] = ID;

		}
	}


	std::array<double,6> get_LC_Pos(int start, int end)
	{
		std::array<double,6> fullLoopInfo;
		bool findStartVertexInMatrix=false,findEndVertexInMatrix=false;
		//get the position of start and end point std::vector<std::array<double,4>> VertexInf
		for(std::vector<std::array<double,4>>::const_iterator itVertex = VertexInf.begin(), 
			lendVertex = VertexInf.end();itVertex!=lendVertex;itVertex++)
		{
			if(int((*itVertex)[0])==start)
			{
                   fullLoopInfo[1]=(*itVertex)[1]; fullLoopInfo[2]=(*itVertex)[2]; //fullLoopInfo[2]=(*itVertex)[3];
				findStartVertexInMatrix=true;

			}
			if(int((*itVertex)[0])==end)
			{
				fullLoopInfo[4]=(*itVertex)[1]; fullLoopInfo[5]=(*itVertex)[2];// fullLoopInfo[2]=(*itVertex)[3];
				findEndVertexInMatrix=true;
			}
			if(findStartVertexInMatrix and findEndVertexInMatrix)
				break;
		}
		if(!(findStartVertexInMatrix and findEndVertexInMatrix))
			cout<<"can't find the position of the poses of the loop"<<endl;
		if(!findStartVertexInMatrix)
		{
			cout<<"start vertex id: "<<start<<endl;
		}
		if(!findEndVertexInMatrix)
		{
			cout<<"end vertex id: "<<end<<endl;
		}

		fullLoopInfo[0]=start;//  fullLoopInfo[1]=startPosition[0];  fullLoopInfo[2]=startPosition[1];
		fullLoopInfo[3]=end;  //  fullLoopInfo[4]=endPosition[0];    fullLoopInfo[5]=endPosition[1];
		return fullLoopInfo;
	}

	std::array<double,3> find_nearest_cluster(std::array<double,6> fullLoopInfo,  std::vector<cluster>& _clustersFound)
	{
		double dis = 0, ID = 0;
		double id_of_nearestLC = 0;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;

		// for(size_t i=0; i<_clustersFound.size(); i++)
		// // for(size_t i=i<_clustersFound.size()-1; i<_clustersFound.size(); i++)
		// {
		// 	returnmid=_clustersFound[i].cal_distance(fullLoopInfo);
		// 	if(i == 0)
		// 	{
		// 		ID = 0;
		// 		dis = returnmid[1];
		// 		id_of_nearestLC = returnmid[0];
		// 	}
			
		// 	else if(returnmid[1] < dis)
		// 	{
		// 		ID = i;
		// 		dis = returnmid[1];
		// 		id_of_nearestLC = returnmid[0];
		// 	}
		// 	//debug: find cluster 9 error

		// 	if(_clustersFound[i].positionserial.size() - _clustersFound[i].dis_cluster_start_end.size() != 1)
		// 	{
		// 		printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
		// 		exit(0);
		// 	}

		// }
		// returndis[0] = ID;
		// returndis[1] = dis;
		// returndis[2] = id_of_nearestLC;
		// return returndis;

		for(size_t i=i<_clustersFound.size()-1; i<_clustersFound.size(); i++)
		{
			returnmid=_clustersFound[i].cal_distance(fullLoopInfo);

				ID = i;
				dis = returnmid[1];
				id_of_nearestLC = returnmid[0];

			//debug: find cluster 9 error

			if(_clustersFound[i].positionserial.size() - _clustersFound[i].dis_cluster_start_end.size() != 1)
			{
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}

		}
		returndis[0] = ID;
		returndis[1] = dis;
		returndis[2] = id_of_nearestLC;
		return returndis;
	}
		
	std::array<double,2> chi2_test(const std::vector<double> & dis_clu, const std::vector<double> & dis_clu_real)
	{
		std::array<double,2> p_value = {0, 0};
		double dis_backup;

		double sum = std::accumulate(std::begin(dis_clu), std::end(dis_clu), 0.0);  
		double mean =  sum / dis_clu.size(); //均值  
						  
		double accum  = 0.0;  
		std::for_each (std::begin(dis_clu), std::end(dis_clu), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
		double chi_statis = (accum/mean); //卡方统计量
		p_value[0] = 1-chi2_p(dis_clu.size()-1, chi_statis);	
		if(dis_clu.size()>2)
		{
			if(dis_clu_real.size() != dis_clu.size()-1)
			{
				cout<<"dis_clu_real.size: "<<dis_clu_real.size()<<endl;
				cout<<"dis_clu.size: "<<dis_clu.size()<<endl;
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			
			sum = std::accumulate(std::begin(dis_clu_real), std::end(dis_clu_real), 0.0);  
			mean =  sum / dis_clu_real.size(); //均值  
			accum  = 0.0;  
			std::for_each (std::begin(dis_clu_real), std::end(dis_clu_real), [&](const double d) {  
			accum  += (d-mean)*(d-mean);  
			});  
						  
			 chi_statis = (accum/mean); //卡方统计量
			p_value[1] = 1-chi2_p(dis_clu_real.size()-1, chi_statis);	
		}
		else
		{
			p_value[1] = p_value[0];
		}
		// cout<<"p_value:"<<p_value[0]<<"  pre_p_value:"<<p_value[1]<<"  dof:"<<(dis_clu.size()-1)<<"  chi_statis:"<<
		// 	chi_statis<<endl;
		return p_value;
	}

	static bool cmp( std::pair<double, int> p,  std::pair<double, int> q)
	{
		return p.first < q.first;
	}



	void clusterize_zihao( const IntPairSet& loops, const char* filename, double thres)//, std::vector<int>& membership, int& clusterCount)
	{
		double disSTART, disEND,  pre_p_value = 0.5;
		int num_loop=0, check = 0;
		bool findStartVertexInMatrix=false,findEndVertexInMatrix=false, ready2chiTest=false;
		std::array<double,3> startPosition, endPosition, nearest_cluster;
		std::array<double,2> tem_dis, p_value;
		std::array<double,6> fullLoopInfo;//six elements:start ID,X,Y,end ID,X,Y		

		collect_vertex(filename);

		if(loops.empty())
		{
			std::cerr<<"clusterize(): "<<__LINE__<<" no loops to make clusters"<<std::endl;
			return;
		}
		_clustersFound.clear();
		for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		{
			int start 	= std::max(it->first,it->second);
			int end 	= std::min(it->first,it->second);



			//print loop number and cluster id of the loop
			num_loop++;
			// cout<<"loop "<<num_loop<<" "<<start<<" "<<end<<endl;

			fullLoopInfo = get_LC_Pos( start,  end);//get loop closure vextexes position


			if(_clustersFound.empty())
			{
				cluster s(start,end,_clustersFound.size());
				s.positionserial.push_back(fullLoopInfo);

				_clustersFound.push_back(s);

			// cout<<"fullLoopInfo: "<<fullLoopInfo[0]<<" "<<fullLoopInfo[1]<<" "<<fullLoopInfo[2]<<" "<<fullLoopInfo[3]<<
			// " "<<fullLoopInfo[4]<<" "<<fullLoopInfo[5]<<endl;

						cout<<"fullLoopInfo: "<<s.positionserial[0][0]<<" "<<s.positionserial[0][1]<<" "
						<<s.positionserial[0][2]<<" "<<s.positionserial[0][3]<<
			" "<<s.positionserial[0][4]<<" "<<s.positionserial[0][5]<<endl;


				clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
				loopToClusterIDMap[*it] = _clustersFound.size()-1;
			}
			else
			{
				// cluster* currentCluster = NULL;

				//search for the nearest cluster to the loop
				nearest_cluster = find_nearest_cluster(fullLoopInfo, _clustersFound);
				

							// if(fullLoopInfo[0] == 3521)
							// {
							// 	for (int p=0;p<_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size();p++)
							// 	{
							// 		cout<<_clustersFound[nearest_cluster[0]].dis_cluster_start_end[p]<<" "<<endl;
							// 	}
							// 	cout<<"nearestCLusterID "<<nearest_cluster[0]<<" pos id:"<<nearest_cluster[2]<<" distance: "<<nearest_cluster[1]<<
							// 		" num of clusters: "<<_clustersFound.size()<<endl;
							// 		 exit(0);
							// }


				// cout<<"zihao_cluster size of dis_new_LC:"<<_clustersFound[nearest_cluster[0]].dis_new_LC.size()<<endl;
				// cout<<"zihao_cluster size of _clustersFound:"<<_clustersFound.size()<<endl;
				// for(int i=0; i<_clustersFound.size();i++)
				// {
				// 	cout<<"cluster "<<i<<" has "<<_clustersFound[i].dis_new_LC.size()<<"elements in dis_new_LC"<<endl;
				// }

				
				//update distance serial of the nearest cluster
				int inspo = _clustersFound[nearest_cluster[0]].update_distance(round(nearest_cluster[2]),nearest_cluster[1]);
				// cout<<"inspo:"<<inspo<<endl;

				//if the size of dis_cluster is bigger than one, we do chi2 test    dis_cluster_backup
				// if(_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size()==1)
				if(_clustersFound[nearest_cluster[0]].dis_new_LC.size()==1)
				{
					if(_clustersFound[nearest_cluster[0]].positionserial.size()!=1)
					{
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}

					clusterIDtoLoopsMap[nearest_cluster[0]].insert(*it);
					loopToClusterIDMap[*it] = nearest_cluster[0];

					_clustersFound[nearest_cluster[0]].dis_cluster_start_end.push_back(nearest_cluster[1]);
					_clustersFound[nearest_cluster[0]].positionserial.push_back(fullLoopInfo);
					//debug: find out why one loop should have distance element in dis_cluster_start_end
					if(_clustersFound[nearest_cluster[0]].positionserial.size() == 1)
					{
						printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
						cout<<"should not add distance to :"<<_clustersFound[nearest_cluster[0]].id_cluster<<endl;
						exit(0);
					}
					
				}
				else
				{
					// _clustersFound[nearest_cluster[0]].update_distance();	 
					// p_value = chi2_test(_clustersFound[nearest_cluster[0]].dis_cluster_start_end);
					p_value = chi2_test(_clustersFound[nearest_cluster[0]].dis_cluster_backup, 
						_clustersFound[nearest_cluster[0]].dis_cluster_start_end);


					 // cout<<"thres: "<<thres<<endl;

					// if(nearest_cluster[1] < 2 )
					// {
					//  	clusterIDtoLoopsMap[nearest_cluster[0]].insert(*it);
					// 	loopToClusterIDMap[*it] = nearest_cluster[0];
					// 	_clustersFound[nearest_cluster[0]].updateAfterChi2OK(fullLoopInfo, inspo);
					// } //else 
					if((p_value[0] > 0.05) and (abs(p_value[0] - p_value[1]) < thres))//
					 {
					 	
					 	clusterIDtoLoopsMap[nearest_cluster[0]].insert(*it);
						loopToClusterIDMap[*it] = nearest_cluster[0];
						_clustersFound[nearest_cluster[0]].updateAfterChi2OK(fullLoopInfo, inspo);
		
					 }
					 else
					 {
					 	int num = _clustersFound[nearest_cluster[0]].dis_cluster_backup.size();
					 	double dis1 = *(_clustersFound[nearest_cluster[0]].dis_new_LC.begin());
					 	double dis2 = *(_clustersFound[nearest_cluster[0]].dis_new_LC.begin()+1);


					 	if(num>2)
					 	{
					 		std::array<double,6> lpb,lpa,lpme;
					 		std::vector <double> dis_to_extract, dis_to_extractx;
					 		double disba,disx;
					 		bool cona,conb;

					 		dis_to_extract.assign(_clustersFound[nearest_cluster[0]].dis_cluster_start_end.begin(), _clustersFound[nearest_cluster[0]].dis_cluster_start_end.end());  
					 		dis_to_extractx.assign(dis_to_extract.begin(),dis_to_extract.end());

					 		//yao kao lv zai mo wei de qing juang
					 		lpme = _clustersFound[nearest_cluster[0]].positionserial[nearest_cluster[2]];

					 		if(nearest_cluster[2] > _clustersFound[nearest_cluster[0]].dis_cluster_start_end.size())
					 		{
					 			cout<<"nearest_cluster[2] "<<nearest_cluster[2]<<" _clustersFound.size"<<_clustersFound.size()<<endl;
						 		cout<<_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size()<<endl;
					 			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
								exit(0);
					 		}

					 		if(nearest_cluster[2] == 0)
					 		{
						 		// lpb = _clustersFound[nearest_cluster[0]].positionserial[nearest_cluster[2]-1];
						 		lpa = _clustersFound[nearest_cluster[0]].positionserial[nearest_cluster[2]+1];
					 			cona = _clustersFound[nearest_cluster[0]].dis_cluster_start_end[nearest_cluster[2]] > nearest_cluster[1];
					 			conb = 1;

						 		dis_to_extractx.erase(dis_to_extractx.begin()+nearest_cluster[2]);
						 		// dis_to_extractx.erase(dis_to_extractx.begin()+nearest_cluster[2]-1);
						 		// dis_to_extractx.insert(dis_to_extractx.begin()+nearest_cluster[2]-1,disba);

					 		}
					 		else if(nearest_cluster[2] == (_clustersFound[nearest_cluster[0]].positionserial.size() - 1))
					 		{
						 		lpb = _clustersFound[nearest_cluster[0]].positionserial[nearest_cluster[2]-1];
						 		// lpa = _clustersFound[nearest_cluster[0]].positionserial[nearest_cluster[2]+1];
					 			conb = _clustersFound[nearest_cluster[0]].dis_cluster_start_end.back() > nearest_cluster[1];
					 			cona = 1;

					 			dis_to_extractx.pop_back();
						 		// dis_to_extractx.erase(dis_to_extractx.begin()+nearest_cluster[2]-1);
						 		// dis_to_extractx.insert(dis_to_extractx.begin()+nearest_cluster[2]-1,disba);
					 		} 
					 		else
					 		{
					 			lpb = _clustersFound[nearest_cluster[0]].positionserial[nearest_cluster[2]-1];
					 			lpa = _clustersFound[nearest_cluster[0]].positionserial[nearest_cluster[2]+1];
					 			disba = sqrt(pow(lpb[1]-lpa[1], 2) + pow(lpb[2]-lpa[2], 2)) + sqrt(pow(lpb[4]-lpa[4], 2) + pow(lpb[5]-lpa[5], 2));

					 			cona = _clustersFound[nearest_cluster[0]].dis_cluster_start_end[nearest_cluster[2]] > nearest_cluster[1];
					 			conb = _clustersFound[nearest_cluster[0]].dis_cluster_start_end[nearest_cluster[2]-1] > nearest_cluster[1];

					 			dis_to_extractx.erase(dis_to_extractx.begin()+nearest_cluster[2]);
						 		dis_to_extractx.erase(dis_to_extractx.begin()+nearest_cluster[2]-1);
						 		dis_to_extractx.insert(dis_to_extractx.begin()+nearest_cluster[2]-1,disba);

					 		}
		 		
					 		if(((sqrt(pow(lpme[1]-fullLoopInfo[1], 2) + pow(lpme[2]-fullLoopInfo[2], 2)) + 
					 			sqrt(pow(lpme[4]-fullLoopInfo[4], 2) + pow(lpme[5]-fullLoopInfo[5], 2))) - nearest_cluster[1]) > 0.001)
					 		{
					 			printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
								exit(0);
					 		}
					 		
					 		p_value = chi2_test(dis_to_extract, dis_to_extractx);
					 		if(cona and conb and (p_value[1] > 0.05))
					 		{
					 			cout<<"pvalue[0]: "<<p_value[0]<<" pvalue[1]: "<<p_value[1]<<endl;
								//create new cluster
						 		cluster s(start,end,_clustersFound.size());
						 		//add distance to new cluster
								s.dis_cluster_start_end.push_back(nearest_cluster[1]);
								//add posiinfo to new cluster
								s.positionserial.push_back(*(_clustersFound[nearest_cluster[0]].positionserial.begin()+nearest_cluster[2]));
								s.positionserial.push_back(fullLoopInfo);

								//add new cluster to _clustersFound
								_clustersFound.push_back(s);
								//set LC and ID map
								clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
								loopToClusterIDMap[*it] = _clustersFound.size()-1;
								//handle the ID and LC map of another loop 
								std::pair<int,int> tochange;
								int startID = lpme[0];
								int endID = lpme[3];
								tochange.first = startID;
								tochange.second = endID;
								//delete the last two elements in positionserial and two all elements in dis_serial
								if(nearest_cluster[2] >= _clustersFound[nearest_cluster[0]].positionserial.size())
								{
					 				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
									exit(0);
								}
								if(1 == _clustersFound[nearest_cluster[0]].positionserial.size())
								{
					 				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
									exit(0);
								}
								
								if(nearest_cluster[2] == (_clustersFound[nearest_cluster[0]].positionserial.size() - 1))
								{
									printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
									_clustersFound[nearest_cluster[0]].dis_cluster_start_end.erase(
										_clustersFound[nearest_cluster[0]].dis_cluster_start_end.begin()+nearest_cluster[2]-1);
								}
								else
								{
	
									cout<<"nearest_cluster[2]: "<<nearest_cluster[2]<<endl;
									cout<<" _clustersFound[nearest_cluster[0]].dis_cluster_start_end.size: "<<
										_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size()<<endl;

									cout<<" _clustersFound[nearest_cluster[0]].positionserial.size: "<<
										_clustersFound[nearest_cluster[0]].positionserial.size()<<endl;

									_clustersFound[nearest_cluster[0]].dis_cluster_start_end.erase(
										_clustersFound[nearest_cluster[0]].dis_cluster_start_end.begin()+nearest_cluster[2]);


								}
	
								_clustersFound[nearest_cluster[0]].positionserial.erase(_clustersFound[nearest_cluster[0]].positionserial.begin()+nearest_cluster[2]);


								cout<<"fullloop[0]: "<<fullLoopInfo[0]<<endl;

								if(nearest_cluster[2] > 0 and nearest_cluster[2] < (_clustersFound[nearest_cluster[0]].positionserial.size() - 1))
									_clustersFound[nearest_cluster[0]].dis_cluster_start_end[nearest_cluster[2]-1] = disba;


								if(clusterIDtoLoopsMap[nearest_cluster[0]].count(tochange))
								{
									clusterIDtoLoopsMap[nearest_cluster[0]].erase(tochange);
									clusterIDtoLoopsMap[_clustersFound.size()-1].insert(tochange);
									//need to chage: delete map from loop to cluster ID
									loopToClusterIDMap[tochange] = _clustersFound.size()-1;
									cout<<"1 loop "<< tochange.first<<" and "<<tochange.second<<"to cluster id map: "<<
										loopToClusterIDMap[tochange]<<endl;		
								}
								else
								{
									tochange.first = endID;
									tochange.second = startID;
									
									cout<<"2 else  loop "<< tochange.first<<" and "<<tochange.second<<"to cluster id map: "<<
										loopToClusterIDMap[tochange]<<endl;
		
								    for(std::set<std::pair<int,int >>::iterator it = clusterIDtoLoopsMap[0].begin();
								    	it!=clusterIDtoLoopsMap[0].end();it++)  
								    {  
								            cout << (*it).first <<" "<< (*it).second << "}\n";  
								    } 

									if(!clusterIDtoLoopsMap[nearest_cluster[0]].count(tochange))
									{
										cout<<"can't find the previous loop in mapset"<<endl;
										cout<<"id:"<<startID<<" "<<endID<<endl;
										cout<<clusterIDtoLoopsMap[nearest_cluster[0]].size()<<endl;
										cout<<"total cluster size:"<<_clustersFound.size()<<endl;
										exit(0);
									}
									clusterIDtoLoopsMap[nearest_cluster[0]].erase(tochange);
									clusterIDtoLoopsMap[_clustersFound.size()-1].insert(tochange);
									loopToClusterIDMap[tochange] = _clustersFound.size()-1;

								}

					 		}
					 		else
					 		{
						 		cluster s(start,end,_clustersFound.size());
								s.positionserial.push_back(fullLoopInfo);

								_clustersFound.push_back(s);
								clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
								loopToClusterIDMap[*it] = _clustersFound.size()-1;
					 		}

					 	}
					 	else if((num == 2) and (_clustersFound[nearest_cluster[0]].dis_cluster_start_end[0]) <= nearest_cluster[1])
					 	{
					 		cluster s(start,end,_clustersFound.size());
							s.positionserial.push_back(fullLoopInfo);

							_clustersFound.push_back(s);
							clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
							loopToClusterIDMap[*it] = _clustersFound.size()-1;

					 	}
					 	else 
					 	{
					 		int addPointer = 0;
					 		if(dis2 >= dis1)
					 			addPointer = 0;
					 		else 
					 			addPointer = 1;
					 		//create new cluster
					 		cluster s(start,end,_clustersFound.size());
					 		//add distance to new cluster
	
							s.dis_cluster_start_end.push_back(min(dis1,dis2));
							//add posiinfo to new cluster
							s.positionserial.push_back(*(_clustersFound[nearest_cluster[0]].positionserial.begin()+addPointer));
							s.positionserial.push_back(fullLoopInfo);

							//add new cluster to _clustersFound
							_clustersFound.push_back(s);
							// //delete old position infor and dis info from old cluster
							// cout<<"_clustersFound[nearest_cluster[0]].positionserial.size: "<<_clustersFound[nearest_cluster[0]].positionserial.size()<<endl;
							// cout<<"_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size: "<<
							// 	_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size()<<endl;
							//set LC and ID map
							clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
							loopToClusterIDMap[*it] = _clustersFound.size()-1;
							//handle the ID and LC map of another loop 
							std::pair<int,int> tochange;
							int startID = (*(_clustersFound[nearest_cluster[0]].positionserial.begin()+addPointer))[0];
							int endID = (*(_clustersFound[nearest_cluster[0]].positionserial.begin()+addPointer))[3];
							tochange.first = startID;
							tochange.second = endID;
							//delete the last two elements in positionserial and two all elements in dis_serial
							_clustersFound[nearest_cluster[0]].positionserial.erase(_clustersFound[nearest_cluster[0]].positionserial.begin()+addPointer);
							_clustersFound[nearest_cluster[0]].dis_cluster_start_end.pop_back();
							if(_clustersFound[nearest_cluster[0]].dis_cluster_start_end.size() != 0 or 
								_clustersFound[nearest_cluster[0]].positionserial.size() != 1)
							{
								printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
								exit(0);
							}							//debug
							// std::array<double,6> mammamlai;
							// mammamlai = (*(_clustersFound[nearest_cluster[0]].positionserial.end()-2));
							// cout<<mammamlai[0]<<" "<<mammamlai[1]<<" "<<mammamlai[2]<<endl;
							// cout<<"nearest_cluster:"<<nearest_cluster[0]<<" "<<nearest_cluster[1]<<endl;
							// cout<<"higher id:"<<startID<<" "<<endID<<endl;

							if(clusterIDtoLoopsMap[nearest_cluster[0]].count(tochange))
							{
								clusterIDtoLoopsMap[nearest_cluster[0]].erase(tochange);
								clusterIDtoLoopsMap[_clustersFound.size()-1].insert(tochange);
								//need to chage: delete map from loop to cluster ID
								loopToClusterIDMap[tochange] = _clustersFound.size()-1;
								cout<<"1 loop "<< tochange.first<<" and "<<tochange.second<<"to cluster id map: "<<
									loopToClusterIDMap[tochange]<<endl;		
							}
							else
							{
								tochange.first = endID;
								tochange.second = startID;
								
								cout<<"2 else  loop "<< tochange.first<<" and "<<tochange.second<<"to cluster id map: "<<
									loopToClusterIDMap[tochange]<<endl;
	
							    for(std::set<std::pair<int,int >>::iterator it = clusterIDtoLoopsMap[0].begin();
							    	it!=clusterIDtoLoopsMap[0].end();it++)  
							    {  
							            cout << (*it).first <<" "<< (*it).second << "}\n";  
							    } 

								if(!clusterIDtoLoopsMap[nearest_cluster[0]].count(tochange))
								{
									cout<<"can't find the previous loop in mapset"<<endl;
									cout<<"id:"<<startID<<" "<<endID<<endl;
									cout<<clusterIDtoLoopsMap[nearest_cluster[0]].size()<<endl;
									cout<<"total cluster size:"<<_clustersFound.size()<<endl;
									exit(0);
								}
								clusterIDtoLoopsMap[nearest_cluster[0]].erase(tochange);
								clusterIDtoLoopsMap[_clustersFound.size()-1].insert(tochange);
								loopToClusterIDMap[tochange] = _clustersFound.size()-1;

							}
						}
							// _clustersFound[nearest_cluster[0]].dis_cluster_start_end.size();
					 }
					 
				}
								

			}



		}

		//merge cluster
// exit(0);
		ofstream fileStream;  
		

		std::vector<std::vector<int> > consistent_pair_clusterr_real;
		consistent_pair_clusterr_real.clear();
		std::vector<std::vector<std::pair<double, int> > > nearest_dis_sequence;
		detect_cluster_merge(_clustersFound, consistent_pair_clusterr_real, nearest_dis_sequence, thres);
		cout<<"information about cluster to merge: \n"<<endl;
		fileStream.open("pairs of clusters to merge.txt",ios::trunc);
		for(int i = 0; i < consistent_pair_clusterr_real.size(); i++)
		{
			cout<<"cluster["<<i<<"]"<<"  consistent_pair_clusterr_real[i].size(): "<<consistent_pair_clusterr_real[i].size()<<endl;;
		}
		for(int i = 0; i < consistent_pair_clusterr_real.size(); i++)
		{
			cout<<" "<<endl;
			cout<<" "<<endl;
			if(consistent_pair_clusterr_real[i].size() > 0)
			{
								
				int cy;
				fileStream<<i<<"  ";

				cout<<"cluster["<<i<<"]"<<"'s potential merge cluster: "<<" ";
				for(std::vector<int>::const_iterator itVertex = consistent_pair_clusterr_real[i].begin(), 
					lendVertex = consistent_pair_clusterr_real[i].end();itVertex!=lendVertex;itVertex++)
				{
					cy = *itVertex;
					fileStream<<cy<<"  ";
					cout<<cy<<"  ";
				}	
				cout<<endl;
				cout<<" "<<endl;
				cout<<"cluster["<<i<<"]"<<"'s nearest cluster sequence: "<<" ";
				for(std::vector<std::pair<double, int> >::const_iterator itVertex = nearest_dis_sequence[i].begin(), 
					lendVertex = nearest_dis_sequence[i].end();itVertex!=lendVertex;itVertex++)
				{
					cout<<(*itVertex).second<<"  ";
				}	
				cout<<endl;
			}
			fileStream<<"\n";
		}
		
		fileStream.close();
		//store file
		merge_cluster( consistent_pair_clusterr_real);

		// exit(0);

		std::array<double,6> ty={1,1,1,1,1,1};

		// fileStream.open("clusterFile.g2o",ios::trunc);
		ofstream fileStreamr; 
		fileStreamr.open(nameofclusterfile,ios::trunc);

		for(size_t i=0; i< _clustersFound.size(); i++)
		{
			for(std::vector<std::array<double,6>>::const_iterator itVertex = _clustersFound[i].positionserial.begin(), 
				lendVertex = _clustersFound[i].positionserial.end();itVertex!=lendVertex;itVertex++)
			{
				ty = *itVertex;

				fileStreamr<<i<<"  "<<ty[0]<<"  "<<ty[1]<<"  "<<ty[2]<<"  "<<ty[3]<<"  "<<ty[4]<<"  "<<ty[5]<<"\n";
				// fileStream<<trystdarray[0]<<"\n";
			}	
		}
		fileStreamr.close();
		cout<<"origianl cluster file has been saved"<<endl;

		//consistency check
		bool ready_2_intra;
		ready_2_intra = get_cluster_node_scequence( node_sequence, _clustersFound, merge_vector);
		if(ready_2_intra)
		{
			cal_seg_taransform_se2(node_sequence, OdoInf, transSequence_whole);
		}
		intra_loop_pair_consistency_test(node_sequence, transSequence_whole, _clustersFound);



		inter_cons_check(consistent_pair_clusterr_real, transSequence_whole);

		
		exit(0);

		sleep(3);

#if 0
		if(0)
		{
			std::cout<<" \% Clusters formed "<<_clustersFound.size()<<std::endl;
			std::cout<<"limits = [ "<<std::endl;
			for(size_t i=0 ; i< _clustersFound.size() ; i++)
			{
				std::cout<<i<<" -> sz "<<_clustersFound[i].size<<" :: ";
				std::cout<<" "<<_clustersFound[i].startLow<<" "<<_clustersFound[i].startHigh<<" ";
				std::cout<<" "<<_clustersFound[i].endLow<<" "<<_clustersFound[i].endHigh<<std::endl;;

			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;


			std::cout<<"membership =[ ";
			for(size_t i=0; i<membership.size();i++)
			{
				std::cout<<membership[i]<<" ";
			}
			std::cout<<std::endl;
			std::cout<<"]; "<<std::endl;
		}
#endif

	}

// 	/* collect vertex ID and position and angle in vertexInfo
	
// 	 */

	int collect_vertex(const char* filename)
	{

		ifstream fileStream;  

	    string tmp,temp1;  
	    std::array<double,4> verT={0, 0, 0, 0};
	    std::array<double,8> odoedge_element;
	    std::array<double,11> savemid;
	    std::array<double,6> lcedge_element;
	    std::pair<int, int>  lc_vertex_pair;
	    // char* seg;
	    int count = 0,dddlndgsfdgj=0;// 行数计数器  
	    int position = 0, position2 = 0;  
	    double nul;
        if (position != string::npos)  


	    // fileStream.open("B25b_0.500.g2o",ios::in);//ios::in 表示以只读的方式读取文件
	    fileStream.open(filename,ios::in);  
	    if(fileStream.fail())//文件打开失败:返回0  
	    {  
	    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
	    	exit(0);
	        return 0;  
	    }  
	    else//文件存在  
	    {  
		    while(getline(fileStream,tmp,'\n'))//读取一行  
		    {  
		    	position = tmp.find("VERTEX_SE2");
		    	position2 = tmp.find("EDGE_SE2");
		    	// cout<<tmp<<endl; 	
		    	// exit(0);
		
		    	istringstream stringin(tmp);
		    	if (tmp.size() > 0 )
		    	{
					if(position!=string::npos)  
			    	{			    		
			    		for (dddlndgsfdgj=0;dddlndgsfdgj<5;dddlndgsfdgj++)
			    		{
			    			switch(dddlndgsfdgj)
			    			{
		                    	case 0:
			    					stringin >> temp1;
			    					// cout<<temp1<<endl;
			    					break;
			    				case 1:
			    				case 2:
			    				case 3:
			    				case 4:
			    					stringin >> verT[dddlndgsfdgj-1];
			    					break;
			    				default:
			    					break;
			    			}
			    		}
			    		count++; 
			    		VertexInf.push_back(verT);
			    		// cout<<verT[0]<<" "<<verT[1]<<" "<<verT[2]<<" "<<verT[3]<<endl;
			    		// cout<<VertexInf.size()<<endl;
			    		// cout<<VertexInf[0][0]<<" "<<VertexInf[0][1]<<" "<<VertexInf[0][2]<<" "<<VertexInf[0][3]<<endl;
			    		// exit(0);
			    	}
			    	else if(position2 !=string::npos)
			    	{
			    		bool odobit = 0;

						for (dddlndgsfdgj = 0; dddlndgsfdgj < 12; dddlndgsfdgj++)
			    		{	
			    			switch(dddlndgsfdgj)
			    			{
		                    	case 0:
			    					stringin >> temp1;
			    					// cout<<temp1<<endl;
			    					break;
			    				case 1:
			    					stringin >> lc_vertex_pair.first;
			    					break;
			    				case 2:
			    					stringin >> lc_vertex_pair.second;
			    					if(std::abs(lc_vertex_pair.first - lc_vertex_pair.second) == 1)
			    						odobit = 1;
			    					break;
			    				case 3:
			    					if(odobit)//if its a odometry edge put in vector[8]
			    					{
			    						odoedge_element[0] = lc_vertex_pair.first;
			    						odoedge_element[1] = lc_vertex_pair.second;
			    						stringin >> odoedge_element[2];
			    					}
			    					else
			    						stringin >> lcedge_element[0];
			    					break;
			    				case 4:	
			    				case 5:
			    				case 6:	
			    					if(odobit)//if its a odometry edge put in vector[8]
			    						stringin >> odoedge_element[dddlndgsfdgj-1];
			    					else
			    						stringin >> lcedge_element[dddlndgsfdgj-3];
			    					break;
			    				case 7:
			    				case 8:
			    				case 10:
			    						stringin >> nul;
			    						if(nul != 0)
			    						{
			    							cout<<nul<<" "<<lc_vertex_pair.first<<" "<<lc_vertex_pair.second<<endl;
			    							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
											exit(0);
			    						}
			    					break;		    				
			    				case 9:
			    					if(odobit)//if its a odometry edge put in vector[8]
			    						stringin >> odoedge_element[dddlndgsfdgj-3];
			    					else
			    						stringin >> lcedge_element[dddlndgsfdgj-5];
			    					break;
			    				case 11:
			    					if(odobit)//if its a odometry edge put in vector[8]
			    						stringin >> odoedge_element[dddlndgsfdgj-4];
			    					else
			    						stringin >> lcedge_element[dddlndgsfdgj-6];
			    					break;
			    				default:
			    				{
			    					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
									exit(0);
			    					break;
			    				}
			    			}
			    		}
			    		count++; 
			    		if(odobit)
			    			OdoInf.push_back(odoedge_element);
			    		else
			    		{
			    			LC_Inf[lc_vertex_pair] = lcedge_element;
				    		cout<<lc_vertex_pair.first<<" "<<lc_vertex_pair.second<<" "<<lcedge_element[0]<<" "
				    			<<lcedge_element[1]<<" "<<lcedge_element[2]<<" "<<lcedge_element[3]<<" "<<lcedge_element[4]<<" "
			    				<<lcedge_element[5]<<endl;		
			    			if(lc_vertex_pair.first <lc_vertex_pair.second){
			    				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
								exit(0);
			    			}
			    		}

			    	}
		    	} 
		    } 


		    if(OdoInf.size() != VertexInf.size()-1)
		    {
		    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
		    }
		    for(int i=0;i<OdoInf.size();i++)
		    {
		    	if(std::abs(OdoInf[i][0] - OdoInf[i][1]) != 1 )
		    	{
		    		cout<<"OdoInf[i][0]: "<<OdoInf[i][0]<<" OdoInf[i][1]: "<<OdoInf[i][1]<<endl;
		    		cout<<"i: "<<i<<endl;
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					// exit(0);
		    	}
		    	else if(OdoInf[i][0] != i)
		    	{
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					// exit(0);	
		    	}
		    	// cout<<OdoInf[i][0]<<" "<<OdoInf[i][1]<<" "<<OdoInf[i][2]<<" "<<OdoInf[i][3]<<" "<<OdoInf[i][4]<<" "
		    	// 	<<OdoInf[i][5]<<" "<<OdoInf[i][6]<<" "<<OdoInf[i][7]<<endl;
		    }

		    printf("loop count:%d \n",count);
		    printf("sum:%d \n",int(VertexInf.size()+LC_Inf.size()+OdoInf.size()));//+LC_Inf.size() +OdoInf.size()
		    cout<<"VertexInf.size():"<<VertexInf.size()<<endl;
		    fileStream.close();  
		    // exit(0);
		    return count ;  
		}  
	}


	void cal_seg_taransform_se2(std::vector<std::vector< std::vector<int> > > & node_sequence, std::vector<std::array<double,8>> & OdoInf, 
		std::vector<std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > > &transSequence_whole)
	{

		g2o::SE2 edge1, edge2, edge0;
		g2o::Vector3 mid_vector3;
		std::pair<int, int> seg_pair;
		std::pair<g2o::SE2, Matrix3d> transSeg;
		std::vector<std::pair<g2o::SE2, Matrix3d> > transSequence;
		std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > transSequence_l;
		Matrix3d m = Matrix3d::Identity();  
		bool sequence_first_2_second;
		int from, to;
		std::array<double, 8> odometryedge_1, odometryedge_2;
		std::array<double, 6> loop_edge;


		transSequence_whole.clear();
		for(int iii = 0; iii < node_sequence.size(); iii++)
		{
			vector<int>::iterator iter=find(killed_cluster.begin(),killed_cluster.end(),iii);
			if(iter != killed_cluster.end())
			{
				transSequence_l.clear();
				transSequence_whole.push_back(transSequence_l);
				continue;
			}

			transSequence_l.clear();
			for(int ii = 0; ii<2; ii++)
			{
				transSequence.clear();
				int to_compar = 0;
				if(node_sequence[iii][ii].size() >= 1)
					to_compar = node_sequence[iii][ii].size()-1;
				else
				{

			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
				
				for(int i = 0; i < to_compar; i++)
				{
					seg_pair.first = node_sequence[iii][ii][i];
					seg_pair.second = node_sequence[iii][ii][i+1];
					if(seg_pair.first < seg_pair.second)
					{
						sequence_first_2_second = 1;
						from = seg_pair.first;
						to = seg_pair.second;
					}
					else
					{
						sequence_first_2_second = 0;
						from = seg_pair.second ;
						to = seg_pair.first;
					}
					int odo_size = to - from;
					if(odo_size < 1)// if the adjacent nodes are the same one, set the transfrom matrix to zero and covariance  to zero
					{
						m = Matrix3d::Identity(); 
						mid_vector3[0] = 0;mid_vector3[1] = 0;mid_vector3[2] = 0;
						edge0.fromVector(mid_vector3);
						m(0,0)= 0; m(1,1)= 0; m(2,2) = 0;
						if(sequence_first_2_second == 0)
							transSeg.first = edge0.inverse();
						else
							transSeg.first = edge0;
						transSeg.second = m;
						transSequence.push_back(transSeg);
						// cout<<"to: "<<to<<endl;
						// cout<<"from: "<<from<<endl;
					 //    printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						// exit(0);		
					}
					else if(odo_size == 1)
					{
						m = Matrix3d::Identity(); 
						mid_vector3[0] = OdoInf[from][2];mid_vector3[1] = OdoInf[from][3];mid_vector3[2] = OdoInf[from][4];
						edge0.fromVector(mid_vector3);
						m(0,0)= 1.0/OdoInf[from][5]; m(1,1)= 1.0/OdoInf[from][6]; m(2,2) = 1.0/OdoInf[from][7];
						if(sequence_first_2_second == 0)
							transSeg.first = edge0.inverse();
						else
							transSeg.first = edge0;
						transSeg.second = m;
						transSequence.push_back(transSeg);
					}
					else
					{
						Matrix3d m1 = Matrix3d::Identity(),  m2 = Matrix3d::Identity(), m_m, J1, J2;
						m1(0,0 )= 1.0/OdoInf[from][5]; m1(1,1)= 1.0/OdoInf[from][6]; m1(2,2) = 1.0/OdoInf[from][7];
						mid_vector3[0] = OdoInf[from][2];mid_vector3[1] = OdoInf[from][3];mid_vector3[2] = OdoInf[from][4];
						edge1.fromVector(mid_vector3);

						for(int j=from+1; j<to; j++)
						{
							m2(0,0 ) = 1.0/OdoInf[j][5]; m2(1,1) = 1.0/OdoInf[j][6]; m2(2,2) = 1.0/OdoInf[j][7];
							mid_vector3[0] = OdoInf[j][2]; mid_vector3[1] = OdoInf[j][3]; mid_vector3[2] = OdoInf[j][4];
							edge2.fromVector(mid_vector3);

							Jacobian_4_edge_propagate(edge1, edge2, J1, J2);
							covariance_propagate(m1, m2, J1, J2, m_m);
							m1 = m_m;

							edge1 *= edge2;//update transform
						}
						if(sequence_first_2_second == 0)
							transSeg.first = edge1.inverse();
						else
							transSeg.first = edge1;
						transSeg.second = m1;
						transSequence.push_back(transSeg);// node level
						// cout<<seg_pair.first<<" "<<seg_pair.second<<endl;
						// cout<<m1<<endl;
		     		}
				}

				transSequence_l.push_back(transSequence);	//cluster level	
				if(ii == 1)
				{
					std::pair<int, int> lppair;

					
					if(node_sequence[iii][0].size() != node_sequence[iii][1].size())//start node length should equal to end node length
					{
						    printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
					}
					if(node_sequence[iii][0].size() >= 1)//if the cluster only has one loop
					{
						transSequence.clear();
						if((node_sequence[iii][0].size() == 1) and (transSequence.size() != 0))
						{
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
						}
						// cout<<"transSequence_l.size: "<<transSequence_l.size()<<endl;
						// exit(0);

						for(int lllpc=0; lllpc < node_sequence[iii][0].size(); lllpc++)
						{
							
							lppair.first = node_sequence[iii][0][lllpc];
							lppair.second = node_sequence[iii][1][lllpc];
							
							loop_edge = LC_Inf[lppair];

							Matrix3d m1 = Matrix3d::Identity();
							m1(0,0 )= 1.0/loop_edge[3]; m1(1,1)= 1.0/loop_edge[4]; m1(2,2) = 1.0/loop_edge[5];
							mid_vector3[0] = loop_edge[0];mid_vector3[1] = loop_edge[1];mid_vector3[2] = loop_edge[2];
							edge1.fromVector(mid_vector3);
							transSeg.first = edge1;
							transSeg.second = m1;
							// cout<<"m covariance: "<<m1<<endl;
							// cout<<edge1.toVector()<<endl;

							transSequence.push_back(transSeg);							
						}
						transSequence_l.push_back(transSequence);
					}
					else
					{
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
					}
					// cout<<"iii: "<<iii<<endl;
					// cout<<"transSequence_l[0].size: "<<transSequence_l[0].size()<<endl;
					// cout<<"transSequence_l[1].size: "<<transSequence_l[1].size()<<endl;
					// cout<<"transSequence_l[2].size: "<<transSequence_l[2].size()<<endl;
					// cout<<"node_sequence[iii][0].size(): "<<node_sequence[iii][0].size()<<endl;
					// cout<<"node_sequence[iii][1].size(): "<<node_sequence[iii][1].size()<<endl;
				}
			}
			transSequence_whole.push_back(transSequence_l); //all clusters level
		}
		// exit(0);
	}

	void Jacobian_4_edge_propagate(g2o::SE2 & TransA, g2o::SE2 & TransB, Matrix3d & J1, Matrix3d & J2)
	{
		J1 = Matrix3d::Identity();
		J2 = Matrix3d::Identity();
		J1(0,2) = -TransB[0]*sin(TransA[2])-TransB[1]*cos(TransA[2]);
		J1(1,2) = TransB[0]*cos(TransA[2])-TransB[1]*sin(TransA[2]);
		J2(0,0) = cos(TransA[2]);J2(0,1) = -sin(TransA[2]);J2(1,0) = sin(TransA[2]);J2(1,1) = cos(TransA[2]);
	}

	void covariance_propagate(Matrix3d & cov1, Matrix3d & cov2, Matrix3d & J1, Matrix3d & J2, Matrix3d & result)
	{
		result = J1*cov1*(J1.transpose())+J2*cov2*(J2.transpose());
	}

	std::pair<bool, double> check_single_loop(int loop_to_check, cluster & _clustersFoundi, int loop_checked,
		std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > &transSequence_cluster, int clus_num, int group_num)//, double& statis
	{
		g2o::SE2 loop1, edgeGo, loop2, edgeBack, Edge_midd, Edge_go_midd, Edge_back_midd, transform_interator;
		Matrix3d Cov_loop1, Cov_edgeGo, Cov_loop2, Cov_edgeBack, Cov_midd, Cov_go_midd, Cov_back_midd, J1, J2, Cov_interator;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;

		// std::array<int,8> nodes;
		int go_odo_num = 0, back_odo_num = 0;
		// //nodes loop check
		// nodes[0] = _clustersFoundi.positionserial[loop_checked][0];
		// nodes[1] = _clustersFoundi.positionserial[loop_checked][3];
		if(loop_to_check <= loop_checked)
		{
			cout<<"loop_to_check: "<<loop_to_check<<" loop_checked: "<<loop_checked<<endl;
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		loop1        =  transSequence_cluster[2][loop_checked].first;
		Cov_loop1    =  transSequence_cluster[2][loop_checked].second;
		loop2        =  transSequence_cluster[2][loop_to_check].first;
		Cov_loop2    =  transSequence_cluster[2][loop_to_check].second;
		if(loop_to_check == loop_checked+1)
		{ 
			edgeGo       =  transSequence_cluster[1][loop_checked].first;
			Cov_edgeGo  =  transSequence_cluster[1][loop_checked].second;
			edgeBack     =  transSequence_cluster[0][loop_checked].first;
			Cov_edgeBack =  transSequence_cluster[0][loop_checked].second;	
		}
		else
		{
			edgeGo       =  transSequence_cluster[1][loop_checked].first;
			Cov_edgeGo  =  transSequence_cluster[1][loop_checked].second;
			edgeBack     =  transSequence_cluster[0][loop_checked].first;
			Cov_edgeBack =  transSequence_cluster[0][loop_checked].second;	

			for(int sumodo = 1; sumodo < loop_to_check; sumodo++)
			{

				Matrix3d m2 = Matrix3d::Identity(), m_m, J1_2_odo, J2_2_odo;

				Edge_go_midd  =  transSequence_cluster[1][loop_checked+sumodo].first;
				Cov_go_midd   =  transSequence_cluster[1][loop_checked+sumodo].second;	
				Edge_back_midd  =  transSequence_cluster[0][loop_checked+sumodo].first;
				Cov_back_midd   =  transSequence_cluster[0][loop_checked+sumodo].second;	

				Jacobian_4_edge_propagate(edgeGo, Edge_go_midd, J1, J2);
				covariance_propagate(Cov_edgeGo, Cov_go_midd, J1, J2, m_m);
				Cov_edgeGo = m_m;
				edgeGo *= Edge_go_midd;//update transform

				Jacobian_4_edge_propagate(edgeBack, Edge_back_midd, J1, J2);
				covariance_propagate(Cov_edgeBack, Cov_back_midd, J1, J2, m_m);
				Cov_edgeBack = m_m;
				edgeBack *= Edge_go_midd;//update transform

				if(loop_checked+sumodo+1 == loop_to_check)
					break;
				else if(sumodo == loop_to_check-1)
				{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}
		Jacobian_4_edge_propagate(loop1, edgeGo, J1, J2);//generate jacobian 
		covariance_propagate(Cov_loop1, Cov_edgeGo, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator = loop1 * edgeGo;//update transform

		Edge_midd = loop2.inverse();
		Jacobian_4_edge_propagate(transform_interator, Edge_midd, J1, J2);//generate jacobian 
		covariance_propagate(Cov_interator, Cov_loop2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator *=  Edge_midd;//update transform

		Edge_midd = edgeBack.inverse();
		Jacobian_4_edge_propagate(transform_interator, Edge_midd, J1, J2);//generate jacobian 
		covariance_propagate(Cov_interator, Cov_edgeBack, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
		Cov_interator = Cov_midd;
		transform_interator *=  Edge_midd;//update transform

		// T = transform_interator.toVector();
		Matrix3d mmmm =  Cov_interator.inverse();
		T(0)= transform_interator[0];
		T(1) = transform_interator[1];
		T(2) = transform_interator[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(1) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(2) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);
		cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		cout<<"cov: "<<endl<<Cov_interator<<endl;
		cout<<"transformDistance: "<<transformDistance<<endl;
		cout<<"transform_interator: "<<transform_interator[0]<<" "<<transform_interator[1]<<" "<<transform_interator[2]<<endl;
		cout<<" "<<endl;
		returnV.second = transformDistance;
		if (transformDistance < 7.81)
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;

	}	

		
	void intra_loop_pair_consistency_test(std::vector<std::vector< std::vector<int> > >  & node_sequence,
		std::vector<std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> > > > & transSequence_whole, 
		std::vector<cluster> & _clustersFound)
	{	

		std::pair<int,double> ele_intra_pair;
		std::vector <std::pair<int,double> > ele_intra_group;
		int increase_uncer = 0;
		std::pair<int,int > uncert_pair;
		std::vector<int > inc;
		for(int i = 0; i < _clustersFound.size(); i++)
		{
			increase_uncer = 0;
			_clustersFound[i].uncert_pair_group.clear();
			//basic check
			vector<int>::iterator iter=find(killed_cluster.begin(),killed_cluster.end(),i);
			if(iter != killed_cluster.end())
			{
				inc.push_back(0);
				continue;
			}
			if((_clustersFound[i].positionserial.size() != node_sequence[i][0].size()) )
			{
				cout<<"_clustersFound.size(): "<<_clustersFound.size()<<endl;
				cout<<"node_sequence.size(): "<<node_sequence.size()<<endl;

				cout<<"i: "<<i<<" positionserial.size: "<<_clustersFound[i].positionserial.size()<<" node size: "<<node_sequence[i][0].size()<<endl;
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
			}
			if((node_sequence[i][0].size() != node_sequence[i][1].size()) )
			{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
			}
			if(	((_clustersFound[i].positionserial.size()-1) != transSequence_whole[i][0].size()) )
			{
				cout<<"i: "<<i<<endl;
				cout<<"_clustersFound[i].positionserial.size(): "<<_clustersFound[i].positionserial.size()<<endl;
				cout<<"transSequence_whole[i][0].size: "<<transSequence_whole[i][0].size()<<endl;
				cout<<"transSequence_whole.size: "<<transSequence_whole.size()<<endl;
				cout<<"_clustersFound.size(): "<<_clustersFound.size()<<endl;
				cout<<"node_sequence.size(): "<<node_sequence.size()<<endl;
				cout<<"killed_cluster.size(): "<<killed_cluster.size()<<endl;
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
			}
			if(	(transSequence_whole[i][0].size() != transSequence_whole[i][2].size()-1))
			{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
			}

			_clustersFound[i].consistentGroup.clear();
			ele_intra_group.clear();			

			if(_clustersFound[i].positionserial.size() == 0)
			{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
			}
			for(int secondFor = 0; secondFor < _clustersFound[i].positionserial.size(); secondFor++)
			{

				if(_clustersFound[i].consistentGroup.size() ==0)
				{
					if(secondFor != 0)
					{
						printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);	
					}
					ele_intra_group.clear();
					ele_intra_pair.first = secondFor;
					ele_intra_pair.second = 0;
					ele_intra_group.push_back(ele_intra_pair);
					_clustersFound[i].consistentGroup.push_back(ele_intra_group);
				}
				else 
				{
					int num_cons_group = 0;
					std::vector<std::pair<int, double> > cons_group,uncertain_check;//uncertain_check only used for display when degbug
					std::vector<int > mess;
					std::pair<int,int> check;
					cons_group.clear();
					std::vector<int> processed;
					std::pair<bool, double> consistent2Group;
					for(int thirdFor = 0; thirdFor < _clustersFound[i].consistentGroup.size(); thirdFor++)
					{
						std::vector<int>::iterator iter=find(processed.begin(),processed.end(),_clustersFound[i].consistentGroup[thirdFor].back().first);
						if(iter != processed.end())
						{
							if(_clustersFound[i].consistentGroup.size() < 2)
							{
							printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
							exit(0);
							}	
							cout<<"_clustersFound[i].consistentGroup.size(): "<<_clustersFound[i].consistentGroup.size()<<endl;	
							consistent2Group = check_single_loop(secondFor,_clustersFound[i], 
							_clustersFound[i].consistentGroup[thirdFor][_clustersFound[i].consistentGroup[thirdFor].size()-2].first,transSequence_whole[i], i, thirdFor);
						}
						else
						{
							consistent2Group = check_single_loop(secondFor,_clustersFound[i], 
							_clustersFound[i].consistentGroup[thirdFor].back().first,transSequence_whole[i], i, thirdFor);
							processed.push_back(_clustersFound[i].consistentGroup[thirdFor].back().first);
						}

						if(consistent2Group.first)
						{
							num_cons_group = num_cons_group+1;

							ele_intra_pair.first = thirdFor;
							ele_intra_pair.second = consistent2Group.second;

							cons_group.push_back(ele_intra_pair);	
							uncertain_check.push_back(consistent2Group)	;	
						}
						else if(consistent2Group.second < 20)//11.34
						{

							if(_clustersFound[i].consistentGroup[thirdFor].size() > 1)
							{
								for(int uncer=0; uncer < _clustersFound[i].consistentGroup[thirdFor].size()-1; uncer++)
								{
									consistent2Group = check_single_loop(secondFor,_clustersFound[i], 
										_clustersFound[i].consistentGroup[thirdFor][uncer].first,transSequence_whole[i], i, thirdFor);
									uncertain_check.push_back(consistent2Group)	;
									if(consistent2Group.second < 7.81)
									{
										num_cons_group = num_cons_group+1;

										ele_intra_pair.first = thirdFor;
										ele_intra_pair.second = consistent2Group.second;
										cons_group.push_back(ele_intra_pair);	
										break;
									}
								}
								// cout<<"uncertain check:"<<endl;
								// for(int uncer=0; uncer < uncertain_check.size(); uncer++)
								// {
								// 	cout<<" "<<uncertain_check[uncer].second;
								// }
								// printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
								// exit(0);	
							}
							else
							{
								//label this uncertain one loop, for further verify 
								uncert_pair.first = thirdFor;
								uncert_pair.first = thirdFor+1;
								_clustersFound[i].uncert_pair_group.push_back(uncert_pair);
							}
							
						} 
					}
					if(cons_group.size() != num_cons_group)
					{
						printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						exit(0);
					}
					if(num_cons_group == 0)
					{
						ele_intra_group.clear();
						ele_intra_pair.first = secondFor;
						ele_intra_pair.second = 0;
						ele_intra_group.push_back(ele_intra_pair);
						_clustersFound[i].consistentGroup.push_back(ele_intra_group);
					}
					else if(num_cons_group == 1)
					{
						ele_intra_pair.first = secondFor;
						ele_intra_pair.second = cons_group[0].second;
						_clustersFound[i].consistentGroup[cons_group[0].first].push_back(ele_intra_pair);
					}
					else
					{
						for(int twoCons = 0; twoCons < _clustersFound[i].consistentGroup.size(); twoCons++)
						{
							cout<<"cluster "<<i <<", group "<<twoCons<<endl;
							for(int each_con_group = 0; each_con_group < _clustersFound[i].consistentGroup[twoCons].size(); each_con_group++)
								cout<<_clustersFound[i].consistentGroup[twoCons][each_con_group].first<<" ";
							cout<<endl;
							for(int each_con_group = 0; each_con_group < _clustersFound[i].consistentGroup[twoCons].size(); each_con_group++)
								cout<<_clustersFound[i].consistentGroup[twoCons][each_con_group].second<<" ";
							cout<<endl;	
						}
						cout<<"num_cons_group: "<<num_cons_group<<endl;
								cout<<"uncertain check statis:"<<endl;
								for(int uncer=0; uncer < uncertain_check.size(); uncer++)
								{
									cout<<"group: "<< uncertain_check[uncer].first<<" statis: "<<uncertain_check[uncer].second<<endl;
								}
						// printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						// exit(0);
						for(int assign_multi_cons_group = 0; assign_multi_cons_group < cons_group.size(); assign_multi_cons_group++)
						{
							ele_intra_pair.first = secondFor;
							ele_intra_pair.second = cons_group[assign_multi_cons_group].second;
							_clustersFound[i].consistentGroup[cons_group[assign_multi_cons_group].first].push_back(ele_intra_pair);	
						}
						increase_uncer = increase_uncer + cons_group.size() -1;

					}
				}

			}
			inc.push_back(increase_uncer);
		}
		//check to make sure total numbers of elements in  groups belongs to the same cluster equal to the loop numbers in the cluster 
		if(inc.size() != _clustersFound.size())
		{
			if(inc.size() != _clustersFound.size()-killed_cluster.size())
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}

		}
		for(int i = 0; i < _clustersFound.size(); i++)
		{
			vector<int>::iterator iter=find(killed_cluster.begin(),killed_cluster.end(),i);
			if(iter != killed_cluster.end())
			{
				continue;
			}
			int sumss=0;
			cout<<" "<<endl;
			cout<<"cluster["<<i<<"]"<<" has "<<_clustersFound[i].positionserial.size()<<" loops"<<endl;
			cout<<"which contains "<<_clustersFound[i].consistentGroup.size()<<" consistent groups."<<endl;
			for(int j=0; j < _clustersFound[i].consistentGroup.size(); j++)
			{
				cout<<"group["<<j<<"]"<<" has "<<_clustersFound[i].consistentGroup[j].size()<<" loops"<<endl;
				sumss = sumss+_clustersFound[i].consistentGroup[j].size();
			}
			if(sumss-inc[i] != _clustersFound[i].positionserial.size())
			{
				cout<<"sumss: "<<sumss<<" increase_uncer: "<<" positionserial.size: "<<endl;
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}

	}

	std::pair<bool, double> check_single_loop_inter(std::array<std::pair<g2o::SE2, Matrix3d>, 4 > &transSequence_cluster_inter)//, double& statis
	{
		g2o::SE2 loop1, loop2, Edge_midd;
		Matrix3d Cov1, Cov2, Cov_midd, J1, J2;
		MatrixXd T(1,3),T_inverse(3,1); 
		std::pair<bool, double> returnV;

		loop1        =  transSequence_cluster_inter[0].first;
		Cov1    =  transSequence_cluster_inter[0].second;

		for(int i = 1; i<4; i++)
		{
			loop2   =  transSequence_cluster_inter[i].first;
			Cov2    =  transSequence_cluster_inter[i].second;

			Jacobian_4_edge_propagate(loop1, loop2, J1, J2);//generate jacobian 
			covariance_propagate(Cov1, Cov2, J1, J2, Cov_midd);// update covariance from two covs and two Jacobians
			Cov1 = Cov_midd;
			loop1 = loop1 * loop2;//update transform
		}
		displayCov = Cov1;
		// T = transform_interator.toVector();
		Matrix3d mmmm =  Cov1.inverse();
		T(0)= loop1[0];
		T(1) = loop1[1];
		T(2) = loop1[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(1) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(2) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);

		// double transformDistance = (T *) * T_inverse;
		// double transformDistance = T(0)*T(0)/Cov_interator(0,0) +  T(1)*T(1)/Cov_interator(1,1) +  T(2)*T(2)/Cov_interator(2,2);

		// cout<<"cluster num: "<<clus_num<<" group_num: "<<group_num<<endl<<" loop1: "<<loop_checked<<" loop2: "<<loop_to_check<<endl;
		// cout<<"cov: "<<endl<<Cov_interator<<endl;
		// cout<<"transformDistance: "<<transformDistance<<endl;
		// cout<<"transform_interator: "<<transform_interator[0]<<" "<<transform_interator[1]<<" "<<transform_interator[2]<<endl;
		// cout<<" "<<endl;
		returnV.second = transformDistance;
		if (transformDistance < 7.81)
			returnV.first = 1;
		else
			returnV.first = 0;
		return returnV;

	}

	void inter_cons_check(std::vector<std::vector<int> > & consistent_pair_clusterr_real, 
		std::vector<std::vector<std::vector<std::pair<g2o::SE2, Matrix3d> >  > > &transSequence_cluster)
	{
		std::vector<std::vector<std::pair<std::pair<int,int>, std::pair<int, std::vector<int> > > > > final_cons_group_whole;
		std::vector<std::pair<std::pair<int,int>, std::pair<int, std::vector<int> > > > final_cons_group;
		std::pair<int, std::vector<int> > consistent_in_next_cluster;
		std::pair<int,int> current_culsterID_groupID;
		std::pair<std::pair<int,int>, std::pair<int, std::vector<int> > > ele_of_final;

		std::vector<std::vector<std::pair<int,int> > > nodessequence2inter, nodessequence2inter_start;
		std::pair<int,int> loopNode;
		std::pair<bool, double> retuValue;
		std::vector<std::pair<int,int> > cluster_group_node, cluster_group_node_start;
		std::array<std::pair<g2o::SE2, Matrix3d>, 2 >  transSequence;
		std::array<int,4> nodes;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4 >  transSequence_cluster_inter;

		cout<<"start of segmentation fault"<<endl;
		for(int i = 0; i < merge_vector.size(); i++)
		{
			cout<<" 111 "<<endl;
			cluster_group_node.clear();
			cluster_group_node_start.clear();
			cout<<" i: "<<i<<endl;
			cout<<"consistent_pair_clusterr_real.size: "<<merge_vector.size()<<endl;
			cout<<"consistent_pair_clusterr_real[i].size: "<<merge_vector[i].size()<<endl;
			int clusterID = merge_vector[i][0];//
			cout<<" clusterID: "<<clusterID<<endl;
			for(int j=0; j<_clustersFound[clusterID].consistentGroup.size(); j++)
			{
				cout<<" j: "<<j<<endl;
				int id4loop = _clustersFound[clusterID].consistentGroup[j].back().first;
				cout<<" id4loop: "<<id4loop<<endl;
				loopNode.first = _clustersFound[clusterID].positionserial[id4loop][0];
				loopNode.second = _clustersFound[clusterID].positionserial[id4loop][3];
				cluster_group_node.push_back(loopNode);

				id4loop = _clustersFound[clusterID].consistentGroup[j][0].first;
				loopNode.first = _clustersFound[clusterID].positionserial[id4loop][0];
				loopNode.second = _clustersFound[clusterID].positionserial[id4loop][3];
				cluster_group_node_start.push_back(loopNode);
			}
			nodessequence2inter.push_back(cluster_group_node);
			nodessequence2inter_start.push_back(cluster_group_node);
		}
		// exit(0);

		for(int i = 0; i < nodessequence2inter.size()-1; i++)
		{
			int clusterID = merge_vector[i][0];
			int nextClusterID = merge_vector[i+1][0];

			for(int j=0; j<_clustersFound[clusterID].consistentGroup.size(); j++)//get the information of the current cluster's consistent group's last node
			{
				current_culsterID_groupID.first  = clusterID ;
				current_culsterID_groupID.second = j;
				consistent_in_next_cluster.second.clear();
				nextClusterID = merge_vector[i+1][0];
				consistent_in_next_cluster.first = nextClusterID;
				int loopID = _clustersFound[clusterID].consistentGroup[j].back().first;

				int node_loop1_start = nodessequence2inter[i][j].first;
				int node_loop1_end   = nodessequence2inter[i][j].second;
				int node_loop2_start = 0;
				int node_loop2_end   = 0;
				transSequence_cluster_inter[0].first = transSequence_cluster[clusterID][2][loopID].first;
				transSequence_cluster_inter[0].second = transSequence_cluster[clusterID][2][loopID].second;

				for(int k=0; k <_clustersFound[nextClusterID].consistentGroup.size(); k++)//get the information of the next cluster's consistent group's first node
				{
					node_loop2_start = nodessequence2inter_start[i+1][k].first;
					node_loop2_end   = nodessequence2inter_start[i+1][k].second;
					nodes[0] = node_loop1_start;
					nodes[1] = node_loop1_end;
					nodes[2] = node_loop2_start;
					nodes[3] = node_loop2_end;

					cal_odo_seg(nodes, transSequence);

					int loopID2 = _clustersFound[nextClusterID].consistentGroup[k][0].first;
					transSequence_cluster_inter[1] = transSequence[1];
					transSequence_cluster_inter[2].first  = transSequence_cluster[nextClusterID][2][loopID2].first.inverse();
					transSequence_cluster_inter[2].second = transSequence_cluster[nextClusterID][2][loopID2].second;
					transSequence_cluster_inter[3].first  = transSequence[0].first.inverse();
					transSequence_cluster_inter[3].second = transSequence[0].second;

					retuValue = check_single_loop_inter(transSequence_cluster_inter);

					if(retuValue.first == 1)
					{
						consistent_in_next_cluster.first  = nextClusterID;//the cluster id of the next one
						consistent_in_next_cluster.second.push_back(k);// = k;
						cout<<"displayCov diagnal: "<<endl<<displayCov<<endl;
					}
					else if(retuValue.second <= 30)
					{
				
						if(_clustersFound[clusterID].consistentGroup[j].size() > 1)
						{
							int size    =   _clustersFound[clusterID].consistentGroup[j].size() - 2;
							int id4loop =  _clustersFound[clusterID].consistentGroup[j][size].first;

							nodes[0] = _clustersFound[clusterID].positionserial[id4loop][0];
							nodes[1] = _clustersFound[clusterID].positionserial[id4loop][3];
							cal_odo_seg(nodes, transSequence);

							transSequence_cluster_inter[0] = transSequence_cluster[clusterID][2][id4loop];
							transSequence_cluster_inter[1] = transSequence[1];
							transSequence_cluster_inter[2].first  = transSequence_cluster[nextClusterID][2][loopID2].first.inverse();
							transSequence_cluster_inter[2].second = transSequence_cluster[nextClusterID][2][loopID2].second;
							transSequence_cluster_inter[3].first  = transSequence[0].first.inverse();
							transSequence_cluster_inter[3].second = transSequence[0].second;

							retuValue = check_single_loop_inter(transSequence_cluster_inter);
							if(retuValue.first == 1)
							{
								consistent_in_next_cluster.first  = nextClusterID;//the cluster id of the next one
								consistent_in_next_cluster.second.push_back(k);
								cout<<"displayCov diagnal: "<<endl<<displayCov<<endl;
								int backdd =0;
							}

						}	

					}
					if(k == _clustersFound[nextClusterID].consistentGroup.size()-1)
					{
						if(consistent_in_next_cluster.second.size() == 0)
						{
							int addd = 2;
							while((i+addd) < (nodessequence2inter_start.size()-2))
							{
								// cout<<"i: "<<i<<endl;
								// cout<<"addd: "<<addd<<endl;
								// cout<<"i+addd: "<<i+addd<<endl;

								for(int sek =0; sek < nodessequence2inter_start[i+addd].size(); sek++)
								{
									cout<<"sek: "<<sek<<endl;
									nextClusterID = merge_vector[i+addd][0];
									node_loop2_start = nodessequence2inter_start[i+addd][sek].first;
									node_loop2_end   = nodessequence2inter_start[i+addd][sek].second;
									nodes[0] = node_loop1_start;
									nodes[1] = node_loop1_end;
									nodes[2] = node_loop2_start;
									nodes[3] = node_loop2_end;

									cal_odo_seg(nodes, transSequence);

									int loopID2 = _clustersFound[nextClusterID].consistentGroup[sek][0].first;
									transSequence_cluster_inter[1] = transSequence[1];
									transSequence_cluster_inter[2].first  = transSequence_cluster[nextClusterID][2][loopID2].first.inverse();
									transSequence_cluster_inter[2].second = transSequence_cluster[nextClusterID][2][loopID2].second;
									transSequence_cluster_inter[3].first  = transSequence[0].first.inverse();
									transSequence_cluster_inter[3].second = transSequence[0].second;

									retuValue = check_single_loop_inter(transSequence_cluster_inter);

									if(retuValue.first == 1)
									{
										consistent_in_next_cluster.first  = nextClusterID;//the cluster id of the next one
										consistent_in_next_cluster.second.push_back(sek);// = k;
										addd = nodessequence2inter_start.size();
										cout<<"displayCov diagnal: "<<endl<<displayCov<<endl;

										break;
									}
									else if(retuValue.second <= 30)
									{
								
										if(_clustersFound[clusterID].consistentGroup[j].size() > 1)
										{
											int size    =   _clustersFound[clusterID].consistentGroup[j].size() - 2;
											int id4loop =  _clustersFound[clusterID].consistentGroup[j][size].first;

											nodes[0] = _clustersFound[clusterID].positionserial[id4loop][0];
											nodes[1] = _clustersFound[clusterID].positionserial[id4loop][3];
											cal_odo_seg(nodes, transSequence);

											transSequence_cluster_inter[0] = transSequence_cluster[clusterID][2][id4loop];
											transSequence_cluster_inter[1] = transSequence[1];
											transSequence_cluster_inter[2].first  = transSequence_cluster[nextClusterID][2][loopID2].first.inverse();
											transSequence_cluster_inter[2].second = transSequence_cluster[nextClusterID][2][loopID2].second;
											transSequence_cluster_inter[3].first  = transSequence[0].first.inverse();
											transSequence_cluster_inter[3].second = transSequence[0].second;

											retuValue = check_single_loop_inter(transSequence_cluster_inter);
											if(retuValue.first == 1)
											{
												consistent_in_next_cluster.first  = nextClusterID;//the cluster id of the next one
												consistent_in_next_cluster.second.push_back(sek);
												addd = nodessequence2inter_start.size();
												cout<<"displayCov diagnal: "<<endl<<displayCov<<endl;
												break;
											}

										}	

									}
								}
								addd = addd+1;
							}
							cout<<"if has terminate"<<endl;
						}
						nextClusterID = merge_vector[i+1][0];
					}

				}
				ele_of_final.first = current_culsterID_groupID;
				ele_of_final.second = consistent_in_next_cluster;
				final_cons_group.push_back(ele_of_final);
			}
			
			final_cons_group_whole.push_back(final_cons_group);
			final_cons_group.clear();
		}

		//code to propose the lase cluster

		std::vector<std::pair<int,int> > breakupCluster;
		for(int i = 0; i < final_cons_group_whole.size(); i++)
		{
			for(int j = 0; j < final_cons_group_whole[i].size(); j++)
			{

				if(final_cons_group_whole[i][j].second.second.size() == 0)
					breakupCluster.push_back(final_cons_group_whole[i][j].first );
				cout<<" cluster["<<final_cons_group_whole[i][j].first.first<<"] "<<"group["<<
					final_cons_group_whole[i][j].first.second<<"] has consistent group: ";
				for(int k = 0; k < final_cons_group_whole[i][j].second.second.size(); k++)
					cout<<final_cons_group_whole[i][j].second.second[k]<<" ";
				cout<<" in cluster["<<final_cons_group_whole[i][j].second.first<<"]"<<endl;
				cout<<" "<<endl;
			}
		}

		for(int i = 0; i < breakupCluster.size(); i++)
		{
			cout<<" cluster["<<breakupCluster[i].first<<"] "<<"group["<<
					breakupCluster[i].second<<"] has no  group"<<endl;
		}

		for(int i = 0; i < breakupCluster.size(); i++)
		{
			int cluster_lo = breakupCluster[i].first;
			int group   = breakupCluster[i].second;
			int j;
			for(j = 0; j< 30; j++)
			{
				if(merge_vector[j][0] == cluster_lo)
					break;
				if(j == merge_vector.size()-1)
				{
					printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}
			}
		}

	}




	void cal_odo_seg(std::array<int,4> &nodes, std::array<std::pair<g2o::SE2, Matrix3d>, 2 > &transSequence)
	{
		g2o::SE2 edge1, edge2, edge0;
		g2o::Vector3 mid_vector3;
		std::pair<int, int> seg_pair;
		std::pair<g2o::SE2, Matrix3d> transSeg;
		Matrix3d m = Matrix3d::Identity(),m_m = Matrix3d::Identity();  
		bool sequence_first_2_second;
		int from, to;
		std::array<double, 8> odometryedge_1, odometryedge_2;
		std::array<double, 6> loop_edge;
		int seri = 0;
				
				for(int i = 0; i < 2; i++)
				{
					if(i == 0)
					{
						seg_pair.first = nodes[0];
						seg_pair.second = nodes[2];
					}
					else
					{
						seg_pair.first = nodes[1];
						seg_pair.second = nodes[3];
					}

					if(seg_pair.first < seg_pair.second)
					{
						sequence_first_2_second = 1;
						from = seg_pair.first;
						to = seg_pair.second;
					}
					else
					{
						sequence_first_2_second = 0;
						from = seg_pair.second ;
						to = seg_pair.first;
					}
					int odo_size = to - from;
					if(odo_size < 1)// if the adjacent nodes are the same one, set the transfrom matrix to zero and covariance  to zero
					{
						m = Matrix3d::Identity(); 
						mid_vector3[0] = 0;mid_vector3[1] = 0;mid_vector3[2] = 0;
						edge0.fromVector(mid_vector3);
						m(0,0)= 0; m(1,1)= 0; m(2,2) = 0;
						if(sequence_first_2_second == 0)
							transSeg.first = edge0.inverse();
						else
							transSeg.first = edge0;
						transSeg.second = m;
						transSequence[seri] = transSeg;
						seri++;
						// cout<<"to: "<<to<<endl;
						// cout<<"from: "<<from<<endl;
					 //    printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
						// exit(0);		
					}
					else if(odo_size == 1)
					{
						m = Matrix3d::Identity(); 
						mid_vector3[0] = OdoInf[from][2];mid_vector3[1] = OdoInf[from][3];mid_vector3[2] = OdoInf[from][4];
						edge0.fromVector(mid_vector3);
						m(0,0)= 1.0/OdoInf[from][5]; m(1,1)= 1.0/OdoInf[from][6]; m(2,2) = 1.0/OdoInf[from][7];
						if(sequence_first_2_second == 0)
							transSeg.first = edge0.inverse();
						else
							transSeg.first = edge0;
						transSeg.second = m;
						transSequence[seri] = transSeg;
						seri++;
					}
					else
					{
						Matrix3d m1 = Matrix3d::Identity(),  m2 = Matrix3d::Identity(), m_m, J1, J2;
						m1(0,0 )= 1.0/OdoInf[from][5]; m1(1,1)= 1.0/OdoInf[from][6]; m1(2,2) = 1.0/OdoInf[from][7];
						mid_vector3[0] = OdoInf[from][2];mid_vector3[1] = OdoInf[from][3];mid_vector3[2] = OdoInf[from][4];
						edge1.fromVector(mid_vector3);

						for(int j=from+1; j<to; j++)
						{
							m2(0,0 ) = 1.0/OdoInf[j][5]; m2(1,1) = 1.0/OdoInf[j][6]; m2(2,2) = 1.0/OdoInf[j][7];
							mid_vector3[0] = OdoInf[j][2]; mid_vector3[1] = OdoInf[j][3]; mid_vector3[2] = OdoInf[j][4];
							edge2.fromVector(mid_vector3);

							Jacobian_4_edge_propagate(edge1, edge2, J1, J2);
							covariance_propagate(m1, m2, J1, J2, m_m);
							m1 = m_m;

							edge1 *= edge2;//update transform
						}
						if(sequence_first_2_second == 0)
							transSeg.first = edge1.inverse();
						else
							transSeg.first = edge1;
						transSeg.second = m1;
						transSequence[seri] = transSeg;
						seri++;
						// cout<<seg_pair.first<<" "<<seg_pair.second<<endl;
						// cout<<m1<<endl;
		     		}
				}
	}

	IntPairSet& getClusterByID(int id){
		return clusterIDtoLoopsMap[id];
	}

	size_t clusterCount()
	{
		return _clustersFound.size();
	}

	bool deleteCluster(int clusterID)
	{
		clusterIDtoLoopsMap.erase(clusterID);

		for(IntPairIDMap::iterator it= loopToClusterIDMap.begin();
				it!=loopToClusterIDMap.end(); it++)
		{
			if(it->second == clusterID)
				loopToClusterIDMap.erase(it->first);
		}

		return true;
	}
	void merge_cluster(std::vector<std::vector<int> > & consistent_pair_clusterr_real)
	{
		merge_vector.clear();
		std:vector<int > ele_merge_cluster;
		ele_merge_cluster.clear();
		bool newclu = 0;

		for(int i =0; i < consistent_pair_clusterr_real.size(); i++)
		{
			if(i == consistent_pair_clusterr_real.size()-1)
			{
				if(merge_vector.back().back() == i)
					break;
				else
				{
					ele_merge_cluster.clear();
					ele_merge_cluster.push_back(i);
					merge_vector.push_back(ele_merge_cluster);
					break;
				}
			}
			if(consistent_pair_clusterr_real[i].size() != 0)
			{
				vector<int>::iterator iter=find(consistent_pair_clusterr_real[i].begin(),consistent_pair_clusterr_real[i].end(),i+1);  
					      
				//delete i 
				if(iter!=consistent_pair_clusterr_real[i].end())
				{
					if(consistent_pair_clusterr_real[i+1].size() != 0)
					{
						vector<int>::iterator iter22=find(consistent_pair_clusterr_real[i+1].begin(),consistent_pair_clusterr_real[i+1].end(),i);  
						if(iter22 != consistent_pair_clusterr_real[i+1].end())
							newclu = 0;
						else
							newclu = 1;
					}
					else
						newclu = 1;
				}
				else
					newclu = 1;
			}
			else
				newclu = 1;
			if(newclu == 1)
			{
				ele_merge_cluster.push_back(i);
				merge_vector.push_back(ele_merge_cluster);
				ele_merge_cluster.clear();
			}
			else
				ele_merge_cluster.push_back(i);

		}
		int j = 0;
		for(int i = 0; i < merge_vector.size(); i++)
		{
			j = j + merge_vector[i].size();
			for(int k = 0; k<merge_vector[i].size(); k++)
				cout<<merge_vector[i][k]<<" ";
			cout<<" "<<endl;

		}
		if(j != consistent_pair_clusterr_real.size())
		{
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}

		for(int i = 0; i < merge_vector.size(); i++)
		{
			if(merge_vector[i].size() > 1)
			{
				for(int j = 1; j< merge_vector[i].size(); j++)
				{
					killed_cluster.push_back(merge_vector[i][j]);
					vector<int>::iterator iter=find(killed_cluster.begin(),killed_cluster.end(),i);
		
					_clustersFound[merge_vector[i][0]].positionserial.insert(_clustersFound[merge_vector[i][0]].positionserial.end(),
						_clustersFound[merge_vector[i][j]].positionserial.begin(),_clustersFound[merge_vector[i][j]].positionserial.end());
					_clustersFound[merge_vector[i][j]].positionserial.clear();
		

					clusterIDtoLoopsMap[merge_vector[i][0]].insert(clusterIDtoLoopsMap[merge_vector[i][j]].begin(), clusterIDtoLoopsMap[merge_vector[i][j]].end());

					clusterIDtoLoopsMap.erase(merge_vector[i][j]);

					for(IntPairIDMap::iterator it= loopToClusterIDMap.begin();
						it!=loopToClusterIDMap.end(); it++)
						{
							if(it->second == merge_vector[i][j])
							{
								std::pair<int,int> loop = it->first;
								loopToClusterIDMap.erase(it->first);
								loopToClusterIDMap[loop] = merge_vector[i][0];
							}

						}

				}

			}

			for(int k = 0; k<merge_vector[i].size(); k++)
				cout<<merge_vector[i][k]<<" ";
			cout<<" "<<endl;

		}


	}


};

#endif /* CLUSTER_HPP_ */
