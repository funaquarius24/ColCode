/*********************************************************************************
 *
 * File name:		rec_noc_traffic.h
 * Class name:		RecNOCTraffic
 * Version:         1.6
 *
 * Software:        MCSL Traffic Loader
 * Author:          Zhe Wang (HKUST), Jiang Xu (HKUST), Xiaowen Wu (HKUST), Xuan Wang (HKUST)
 *                  Zhehui Wang (HKUST), Duong Luan (HKUST), Peng Yang (HKUST), 
 *                  Wei Zhang (HKUST), Bin Li (Intel), Ravi Lyer (Intel), Ramesh Illikkal (Intel)
 * Past members:    Weichen Liu, Yaoyao Ye
 * Website:         http://www.ece.ust.hk/~eexu 
 *
 * The copyright information of this program can be found in the file COPYRIGHT.
 *
 *********************************************************************************/

#ifndef MCSL_REC_NOC_TRAFFIC
#define MCSL_REC_NOC_TRAFFIC

#include <assert.h>
#include <string>
#include <iostream>
#include <vector>
#include "rec_proc.h"

using namespace std;

class RecNOCTraffic {

public:
	RecNOCTraffic()		{}
	~RecNOCTraffic()	{}

	// assistant functions
	std::vector<RecProc, std::allocator<RecProc>>					load_traffic(string trace_file);	// load traffic from trace file
	int					print_traffic();					// verify the loaded traffic

	RecProc*			get_proc(int x, int y);				// get processor by coordinate
	RecProc*			get_proc(int id);					// get processor by id
	
	// basic functions
	int					get_topology();
	int					get_num_proc();
	int					get_num_row();
	int					get_num_col();
	int					get_num_task();
	int					get_num_edge();
	int					get_num_iter();

private:
    int                 topology;               	// the topology code
	int					num_row;					// number of rows in mesh/torus
	int					num_col;					// number of columns in mesh/torus
	int                 num_iter;                   // total number of iterations the graph executes for
	vector<RecProc>     proc_list;                  // the list of PBs
	vector<RecTask>     task_list;                  // the list of tasks
	vector<RecEdge>     edge_list;                  // the list of edges

    vector<RecTask*>	starting_task_list;     	// the list of starting tasks
    vector<RecTask*>    finishing_task_list;    	// the list of finishing tasks
};

#endif


