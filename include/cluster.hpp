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

	void find_cons_cluster(IntPairSet::const_iterator & LP_nodes,  std::vector<cluster>& _clustersFound, std::vector<int> & cons_cluster_number,
		std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map, 
		std::vector<std::array<double,4>> & VertexInf, bool debug)
	{
		double dis = 0, ID = 0, id_of_nearestLC = 0, deltaX1, deltaY1, deltaX2, deltaY2, covX, covY;
		int start, end;
		std::array<double,3> returndis;
		std::array<double,2> returnmid;
		std::array<std::pair<g2o::SE2, Matrix3d>, 4> FullInfo;
		g2o::SE2  Trans1, Trans2, midTrans;
		Matrix3d  cov1, cov2, midCov;
		std::pair<int, int> loop_node_pair;
		cons_cluster_number.clear();

		//iterate the elements in clusters, if find one consistent cluster, jump out loop.
		for(int i=_clustersFound.size()-1; i >= 0; i--)
		{
			//get the nearest distance to the cluster to judge whether the consistency check is reliable
			//how to get the vertexinfro 
			if(_clustersFound.size() == 1)
			{
				std::array<double,4> startv = VertexInf[(*LP_nodes).first];
				std::array<double,4> endv = VertexInf[(*LP_nodes).second];

					deltaX1 = abs(_clustersFound[0].positionserial.back()[1] - startv[1]);
					deltaY1 = abs(_clustersFound[0].positionserial.back()[2] - startv[2]);

					deltaX2 = abs(_clustersFound[0].positionserial.back()[4] - endv[1]);
					deltaY2 = abs(_clustersFound[0].positionserial.back()[5] - endv[2]);
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			}
			else
			{
				get_distance_to_neighbor_cluster(i, _clustersFound, deltaX1,  deltaY1, deltaX2,  deltaY2);
				printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
			}
			//synthesize the start odometry edges
			loop_node_pair.first = _clustersFound[i].positionserial.back()[0];
			loop_node_pair.second = _clustersFound[i].positionserial.back()[3];
			start = loop_node_pair.first;
			end   = (*LP_nodes).first;
			if(start > end)
			{
				synthesize_odo_edges( end,  start,  OdoInf, Trans1,  cov1);
				Trans1 = Trans1.inverse();
			}
			else
				synthesize_odo_edges( start,  end,  OdoInf, Trans1,  cov1);
			bool accept = 0;
			if(cov1(0,0) > deltaX1)
				accept = 1;
			if(cov1(1,1) > deltaY1)
				accept = 1;
			//synthesize the end odometry edges
			start = (*LP_nodes).second ;
			end   = loop_node_pair.second;

			if(start > end)
			{
				synthesize_odo_edges( end,  start,  OdoInf, Trans2,  cov2);
				Trans2 = Trans2.inverse();
			}
			else
				synthesize_odo_edges( start,  end,  OdoInf, Trans2,  cov2);
			if(cov2(0,0) > deltaX2)
				accept = 1;
			if(cov2(1,1) > deltaY2)
				accept = 1;
	
			//assign value of 4 parts which from one cycle, the sequence is from part1{lastloop.start > test_loop.start} to part2{test_loop}
			// to part3 {test_loop.second > last_loop.second} to part4{the inverse of last_loop} 
			FullInfo[0] = std::pair<g2o::SE2, Matrix3d> (Trans1, cov1);
			FullInfo[1] = LP_Trans_Covar_Map[(*LP_nodes)];
			FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
			FullInfo[2] = std::pair<g2o::SE2, Matrix3d> (Trans2, cov2);
			midTrans    = LP_Trans_Covar_Map[loop_node_pair].first;
			midTrans = midTrans.inverse();
			midCov      = LP_Trans_Covar_Map[loop_node_pair].second;
			FullInfo[3] = std::pair<g2o::SE2, Matrix3d> (midTrans, midCov);

			std::pair<bool, double> reV;
			double  transX_residual;
			double  transY_residual;
			double  transA_residual;

			reV = check_single_loop_inter(FullInfo, covX, covY, displayCov, transX_residual, transY_residual, transA_residual);
			if(reV.first == 1)
			{
				cout<<"pass chi2"<<endl;
				if(accept == 1)
				{
					cout<<"accept == 1"<<endl;
					cout<<"cov1x: "<<cov1(0,0)<<" cov1Y: "<<cov1(1,1)<<" cov2x: "<<cov2(0,0)<<"cov2Y: "<<cov2(1,1)<<endl;
					cout<<"deltaX1: "<<deltaX1<<" deltaY1: "<<deltaY1<<" deltaX2: "<<deltaX2<<" deltaY2: "<<deltaY2<<endl;
					cout<<"transX_residual: "<<transX_residual<<" transY_residual: "<<transY_residual<<"transA_residual: "<<transA_residual<<endl;
					cout<<"covX: "<<covX<<" covY: "<<covY<<"covAngle: "<<displayCov(2,2)<<" reV.second: "<<reV.second<<endl;
					cout<<"test loop.first: "<<(*LP_nodes).first<<" test loop.second: "<<(*LP_nodes).second<<endl;
					cout<<"loop_node_pair.first: "<<loop_node_pair.first<<" loop_node_pair.second: "<<loop_node_pair.second<<endl;

									cout<<"displayCov:"<<endl;
				cout<<displayCov<<endl;

				cout<<"displayCov.inverse():"<<endl;
				cout<<displayCov.inverse()<<endl;
					cout<<" "<<endl;


					cons_cluster_number.push_back(i);
					continue;
				}
				else
				{
					cout<<"prefect "<<endl;

					cout<<" test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
					cout<<"loop in cluster: "<<loop_node_pair.first<<" "<<loop_node_pair.second<<" in cluster "<< i<<endl;
					cout<<"covX: "<<covX<<" covY: "<<covY<<"covAngle: "<<displayCov(2,2)<<" reV.second: "<<reV.second<<endl;

					cons_cluster_number.push_back(i);
				}
			}
			else
			{
				if(reV.first < 20)
				{

				}

				cout<<" test      loop: "<<(*LP_nodes).first<<" "<<(*LP_nodes).second<<endl;
				cout<<"loop in cluster: "<<loop_node_pair.first<<" "<<loop_node_pair.second<<endl;
				cout<<"cov1x: "<<cov1(0,0)<<" cov1Y: "<<cov1(1,1)<<" cov2x: "<<cov2(0,0)<<"cov2Y: "<<cov2(1,1)<<endl;
				cout<<"deltaX1: "<<deltaX1<<" deltaY1: "<<deltaY1<<" deltaX2: "<<deltaX2<<" deltaY2: "<<deltaY2<<endl;
				cout<<"transX_residual: "<<transX_residual<<" transY_residual: "<<transY_residual<<"transA_residual: "<<transA_residual<<endl;
				cout<<"covX: "<<covX<<" covY: "<<covY<<" covAngle: "<<displayCov(2,2)<<" reV.second: "<<reV.second<<endl;
													cout<<"displayCov:"<<endl;
				cout<<displayCov<<endl;
				cout<<"displayCov.inverse():"<<endl;
				cout<<displayCov.inverse()<<endl;
				cout<<" "<<endl;
			}
		}

	}
	void get_distance_to_neighbor_cluster(int i, std::vector<cluster>& _clustersFound, double & disx1, double & disy1, double & disx2, double & disy2)
	{
		disx1 = 10000;
		disy1 = 10000;
		disx2 = 10000;
		disy2 = 10000;
		for(int j = 0; j<_clustersFound.size(); j++)
		{
			if(j == i)
				continue;
			//disX1
			if(abs(_clustersFound[j].positionserial[0][1] - _clustersFound[i].positionserial.back()[1]) <  disx1)
				disx1 = abs(_clustersFound[j].positionserial[0][1] - _clustersFound[i].positionserial.back()[1]);
			if(abs(_clustersFound[j].positionserial.back()[1] - _clustersFound[i].positionserial.back()[1]) <  disx1)
				disx1 = abs(_clustersFound[j].positionserial.back()[1] - _clustersFound[i].positionserial.back()[1]);

			//disY1
			if(abs(_clustersFound[j].positionserial[0][2] - _clustersFound[i].positionserial.back()[2]) <  disy1)
				disy1 = abs(_clustersFound[j].positionserial[0][2] - _clustersFound[i].positionserial.back()[2]);
			if(abs(_clustersFound[j].positionserial.back()[2] - _clustersFound[i].positionserial.back()[2]) <  disy1)
				disy1 = abs(_clustersFound[j].positionserial.back()[2] - _clustersFound[i].positionserial.back()[2]);
			

			if(abs(_clustersFound[j].positionserial[0][4] - _clustersFound[i].positionserial.back()[4]) <  disx2)
				disx2 = abs(_clustersFound[j].positionserial[0][4] - _clustersFound[i].positionserial.back()[4]);
			if(abs(_clustersFound[j].positionserial.back()[4] - _clustersFound[i].positionserial.back()[4]) <  disx2)
				disx2 = abs(_clustersFound[j].positionserial.back()[4] - _clustersFound[i].positionserial.back()[4]);

			if(abs(_clustersFound[j].positionserial[0][5] - _clustersFound[i].positionserial.back()[5]) <  disy2)
				disy2 = abs(_clustersFound[j].positionserial[0][5] - _clustersFound[i].positionserial.back()[5]);	
			if(abs(_clustersFound[j].positionserial.back()[5] - _clustersFound[i].positionserial.back()[5]) <  disy2)
				disy2 = abs(_clustersFound[j].positionserial.back()[5] - _clustersFound[i].positionserial.back()[5]);
		}
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
		std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	LP_Trans_Covar_Map;	
		//LC_Inf is a map, from nodes pair of loop closure to the six element array of 
		collect_vertexAndLP(filename, LC_Inf, LP_Trans_Covar_Map);

		if(loops.empty())
		{
			std::cerr<<"clusterize(): "<<__LINE__<<" no loops to make clusters"<<std::endl;
			return;
		}
		_clustersFound.clear();

		// //debug check the result of one loop closure multiply the inverse of itself 
		// for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		// {

		// }
		for(IntPairSet::const_iterator it = loops.begin(), lend = loops.end();it!=lend;it++)
		{

			int start 	= it->first;
			int end 	= it->second;
			if(start<end)
			{
				printf("This error about start node and end node is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
			//get loop closure vextexes position
			fullLoopInfo = get_LC_Pos( start,  end);
			//print loop number and cluster id of the loop
			num_loop++;
			cout<<"loop "<<num_loop<<" "<<start<<" "<<end<<endl;
			cout<<"current has "<<_clustersFound.size()<<" clusters"<<endl;

			if(_clustersFound.empty())
			{
				cluster s(start,end,_clustersFound.size());
				s.positionserial.push_back(fullLoopInfo);

				_clustersFound.push_back(s);

						cout<<"fullLoopInfo: "<<s.positionserial[0][0]<<" "<<s.positionserial[0][1]<<" "
						<<s.positionserial[0][2]<<" "<<s.positionserial[0][3]<<
			" "<<s.positionserial[0][4]<<" "<<s.positionserial[0][5]<<endl;


				clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
				loopToClusterIDMap[*it] = _clustersFound.size()-1;
			}
			else
			{
				//search for the nearest cluster to the loop
				std::vector<int>  cons_cluster_number;
		
				find_cons_cluster( it,  _clustersFound,  cons_cluster_number, LP_Trans_Covar_Map, VertexInf, 1);

				//size equal to 0, then it means find no constent cluster, so construct a new cluster
				if(cons_cluster_number.size() == 0)
				{
					cluster s(start,end,_clustersFound.size());
					s.positionserial.push_back(fullLoopInfo);

					_clustersFound.push_back(s);
/*					cout<<"fullLoopInfo: "<<s.positionserial[0][0]<<" "<<s.positionserial[0][1]<<" "
						<<s.positionserial[0][2]<<" "<<s.positionserial[0][3]<<
						" "<<s.positionserial[0][4]<<" "<<s.positionserial[0][5]<<endl;*/
					clusterIDtoLoopsMap[_clustersFound.size()-1].insert(*it);
					loopToClusterIDMap[*it] = _clustersFound.size()-1;
				}
				//find one constent cluster, add the loop to it
				else if (cons_cluster_number.size() == 1)
				{
					cout<<" consistent cluster: "<<cons_cluster_number[0]<<endl;

					// sleep(2);
					int consCluster_serialNum = cons_cluster_number[0];
					_clustersFound[consCluster_serialNum].positionserial.push_back(fullLoopInfo);

					clusterIDtoLoopsMap[consCluster_serialNum].insert(*it);
					loopToClusterIDMap[*it] = consCluster_serialNum;
				}
				//find more than one consistent cluster
				else 
				{
					cout<<"consistent more than one cluster that is consistent to the test loop"<<endl;
					cout<<" cons_cluster_number.size is "<<cons_cluster_number.size()<<endl;
					
					cout<<" cons_cluster_number element:  "<<cons_cluster_number[0]<<" "<<cons_cluster_number[1]<<endl;
					cout<<" numbers in pervious clusters:  "<<_clustersFound[cons_cluster_number[0]].positionserial.size()<<
						" "<<_clustersFound[cons_cluster_number[1]].positionserial.size()<<endl;
					printf("This fake error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);
				}					
			}
		}
		//merge cluster
// exit(0);
		std::array<double,6> ty={1,1,1,1,1,1};

		// fileStream.open("clusterFile.g2o",ios::trunc);
		ofstream fileStreamr; 
		fileStreamr.open(nameofclusterfile,ios::trunc);
		// cout<<"clusters:"<<endl;
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

	int collect_vertexAndLP(const char* filename, std::map<std::pair<int, int>, std::array<double,6> > &LC_Inf,
		std::map<std::pair<int, int>, std::pair<g2o::SE2, Matrix3d> >	& LP_Trans_Covar_Map)
	{

		ifstream fileStream;  

	    string tmp,temp1;  
	    std::array<double,4> verT={0, 0, 0, 0};
	    std::array<double,8> odoedge_element;
	    std::array<double,11> savemid;
	    std::array<double,6> lcedge_element;
	    std::pair<int, int>  lc_vertex_pair;
	    std::pair<g2o::SE2, Matrix3d> ele_lp_trans_covar;
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
			    		
			    		if(odobit)
			    			OdoInf.push_back(odoedge_element);
			    		else
			    		{
			    			count++; 
			    			LC_Inf[lc_vertex_pair] = lcedge_element;
				    		cout<<lc_vertex_pair.first<<" "<<lc_vertex_pair.second<<" "<<lcedge_element[0]<<" "
				    			<<lcedge_element[1]<<" "<<lcedge_element[2]<<" "<<lcedge_element[3]<<" "<<lcedge_element[4]<<" "
			    				<<lcedge_element[5]<<endl;	

							Matrix3d m1 = Matrix3d::Identity();
							g2o::Vector3 mid_vector3; 
							g2o::SE2 edge1;

							m1(0,0 )= 1.0/lcedge_element[3]; 
							m1(1,1)= 1.0/lcedge_element[4]; 
							m1(2,2) = 1.0/lcedge_element[5];

							mid_vector3[0] = lcedge_element[0];
							mid_vector3[1] = lcedge_element[1];
							mid_vector3[2] = lcedge_element[2];
							edge1.fromVector(mid_vector3);

							ele_lp_trans_covar.first  = edge1;
							ele_lp_trans_covar.second = m1;

							LP_Trans_Covar_Map[lc_vertex_pair] = ele_lp_trans_covar;

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
					exit(0);
		    	}
		    	else if(OdoInf[i][0] != i)
		    	{
			    	printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
					exit(0);	
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

	//input: the start and end node serial number correspond to the segment you want to synthesize
	//output: the tranform and covariance info of the synthesized segment, in Trans and cov
	void synthesize_odo_edges(int start, int end, std::vector<std::array<double,8>> & OdoInf, g2o::SE2 & Trans, Matrix3d & cov)
	{
		std::array<double,8> startNodeInfo = OdoInf[start];
		g2o::Vector3  mid_vector3;
		Matrix3d m_m, m2, J1, J2;
		g2o::SE2 edge2;
		int lengthNode = end -start;
		if(startNodeInfo[0] != start){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		if(lengthNode < 0){
			printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
			exit(0);
		}
		cov = Matrix3d::Identity(); 
		if(lengthNode == 0){
			cov(0,0)= 0; cov(1,1)= 0; cov(2,2) = 0;
			mid_vector3[0] = 0;mid_vector3[1] = 0;mid_vector3[2] = 0;
			Trans.fromVector(mid_vector3);
		}
		else
		{
			if(lengthNode == 1 and (startNodeInfo[1] != end))
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}			
			cov(0,0)= 1.0/startNodeInfo[5]; cov(1,1)= 1.0/startNodeInfo[6]; cov(2,2) = 1.0/startNodeInfo[7];
			mid_vector3[0] = startNodeInfo[2];
			mid_vector3[1] = startNodeInfo[3];
			mid_vector3[2] = startNodeInfo[4];
			Trans.fromVector(mid_vector3);
			int j=start;
			for(j=start+1; j<end; j++)
			{
				m2 = Matrix3d::Identity();
				m2(0,0 ) = 1.0/OdoInf[j][5]; m2(1,1) = 1.0/OdoInf[j][6]; m2(2,2) = 1.0/OdoInf[j][7];
				mid_vector3[0] = OdoInf[j][2]; mid_vector3[1] = OdoInf[j][3]; mid_vector3[2] = OdoInf[j][4];
				edge2.fromVector(mid_vector3);

				Jacobian_4_edge_propagate(Trans, edge2, J1, J2);
				covariance_propagate(cov, m2, J1, J2, m_m);
				cov = m_m;

				Trans *= edge2;//update transform
			}
			if(OdoInf[j-1][1] != end)
			{
				printf("This error is in %s on line %d\n",  __FILE__, __LINE__);
				exit(0);
			}
		}
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

//input is four pair<transform, covariance> of two loop and two odo edge segments
//return is one pair<pass_check_or_not, transfrom_distance>
	std::pair<bool, double> check_single_loop_inter(std::array<std::pair<g2o::SE2, Matrix3d>, 4 > &transSequence_cluster_inter, 
	double& covX, double & covY, Matrix3d & displayCov,double & transX_residual, double & transY_residual, double & transA_residual)//, double& statis
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
		covX = Cov1(0, 0);
		covY = Cov1(1, 1);
		// T = transform_interator.toVector();
		Matrix3d mmmm =  Cov1.inverse();
		T(0)= loop1[0];
		T(1) = loop1[1];
		T(2) = loop1[2];	
		T_inverse(0) = T(0) * mmmm(0,0) + T(1) * mmmm(1,0) + T(2) * mmmm(2,0);	
		T_inverse(1) = T(0) * mmmm(0,1) + T(1) * mmmm(1,1) + T(2) * mmmm(2,1);
		T_inverse(2) = T(0) * mmmm(0,2) + T(1) * mmmm(1,2) + T(2) * mmmm(2,2);
		double transformDistance = T_inverse(0)*T(0) + T_inverse(1)*T(1) + T_inverse(2)*T(2);
		transX_residual = T(0);
		transY_residual = T(1);
		transA_residual = T(2);

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
	}
		
};

#endif /* CLUSTER_HPP_ */
