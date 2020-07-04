/*
 * Noxim - the NoC Simulator 
 *
 * (C) 2005-2018 by the University of Catania
 * For the complete list of authors refer to file ../doc/AUTHORS.txt
 * For the license applied to these sources refer to file ../doc/LICENSE.txt
 *
 * This file contains the declaration of the processing element
 */

#ifndef __NOXIMPROCESSINGELEMENT_H__
#define __NOXIMPROCESSINGELEMENT_H__

#include <queue>
#include <systemc.h>

#include "DataStructs.h"
#include "GlobalTrafficTable.h"
#include "Utils.h"
#include "tasks/rec_proc.h"
#include "tasks/rec_task.h"
#include "tasks/rec_edge.h"
#include "tasks/rec_noc_traffic.h"
#include "tasks/stat_proc.h"
#include "tasks/stat_task.h"
#include "tasks/stat_edge.h"
#include "tasks/stat_noc_traffic.h"


using namespace std;

SC_MODULE(ProcessingElement)
{

    // I/O Ports
    sc_in_clk clock;		// The input clock for the PE
    sc_in < bool > reset;	// The reset signal for the PE

    sc_in < Flit > flit_rx;	// The input channel
    sc_in < bool > req_rx;	// The request associated with the input channel
    sc_out < bool > ack_rx;	// The outgoing ack signal associated with the input channel
    sc_out < TBufferFullStatus > buffer_full_status_rx;	

    sc_out < Flit > flit_tx;	// The output channel
    sc_out < bool > req_tx;	// The request associated with the output channel
    sc_in < bool > ack_tx;	// The outgoing ack signal associated with the output channel
    sc_in < TBufferFullStatus > buffer_full_status_tx;

    sc_in < int >free_slots_neighbor;

    // Registers
    int local_id;		// Unique identification number
    bool current_level_rx;	// Current level for Alternating Bit Protocol (ABP)
    bool current_level_tx;	// Current level for Alternating Bit Protocol (ABP)
    queue < Packet > packet_queue;	// Local queue of packets
    bool transmittedAtPreviousCycle;	// Used for distributions with memory

    // Functions
    void rxProcess();		// The receiving process
    void txProcess();		// The transmitting process
    bool canShot(Packet & packet);	// True when the packet must be shot
    Flit nextFlit();	// Take the next flit of the current packet
    Packet trafficTest();	// used for testing traffic
    Packet trafficRandom();	// Random destination distribution
    Packet trafficTranspose1();	// Transpose 1 destination distribution
    Packet trafficTranspose2();	// Transpose 2 destination distribution
    Packet trafficBitReversal();	// Bit-reversal destination distribution
    Packet trafficShuffle();	// Shuffle destination distribution
    Packet trafficButterfly();	// Butterfly destination distribution
    Packet trafficLocal();	// Random with locality
    Packet trafficULocal();	// Random with locality
    void traffic_mcsl();  //MCSL traffic


    // MCSL
    RecProc *rec_proc;
    RecProc local_rec_proc;

    StatProc *stat_proc;
    StatProc local_stat_proc;
    
    //RecNOCTraffic *rec_noc;
    
    RecTask *rec_task;
    StatTask *stat_task;

    vector<RecTask*> *taskList;
    //vector<StatTask*> *taskList;

    // ism Since we broke the whole message into multiple packets, we don't want to
    // send everything at once. so we queue them up before sending
    queue<Packet> packet_gen_queue;
    Packet packet_tmpl;
    int packets_left;
    int packet_number = 0;

    vector<pair<int, int>> task_list_int;
    vector<vector<pair<int, int>>> all_task_list;
    
    vector<TaskInt> list_of_taskInt;
    

    // both incoming dependencies and received edges
    //vector<RecEdge*> *incoming_dep_edge_list;
    //vector<pair<int, int>> incoming_dep_edge_list_int;

    // both outgoing dependencies and received edges
    //vector<RecEdge*> *outgoing_dep_edge_list;
    vector<pair<int, int>> outgoing_dep_edge_list_int;

    /////////////////////////////////////////////////////////
    int incoming_dep_edge_pos;
    int outgoing_dep_edge_pos;
    int task_list_pos;
    ////////////////////////////////////////////////////////

    // ism map that saves all the received edge numbers and their tasks
    map<int, vector<int>> receivedEdgeNumbers;

    // ism can't come up with better name
    // it saves packet's timestamp and how many of them received so far
    map<int, int> map_all_edge_packets;

    // ism map<mcslCount, map<edgeNumber, pair<count, max>>>    
    // where count is the number of edge packets received so far and max is 
    // the max number of packets expected
    map<int, map<int, pair<int, int>>> map_all_edge_packets_itr;

    // Time it takes for task to be executed
    vector<int> *execTimes;
    int execTime;

    RecEdge rec_edge;
    StatEdge stat_edge;
    int taskPointer = -1;
    int edgePointer = -1;
    bool ongoing_task = false;
    bool mcsl_ready = false;
    int recIterationNumber;  // ism Recorded iteration number for this pc
    bool finishedTaskIteration = true;  // ism This allows to determine whether to move to next task
    int startTaskTime = 0;

    
    int mcslCount = -1;  // ism mcslCount is iteration number
    int count; // ism only meant to test for debug purposes...

    bool metDeps();
    void clearDeps(Flit &flit);
    void clearReadymadeDeps();


    GlobalTrafficTable *traffic_table;	// Reference to the Global traffic Table
    bool never_transmit;	// true if the PE does not transmit any packet 
    //  (valid only for the table based traffic)

    void fixRanges(const Coord, Coord &);	// Fix the ranges of the destination
    int randInt(int min, int max);	// Extracts a random integer number between min and max
    int getRandomSize();	// Returns a random size in flits for the packet
    void setBit(int &x, int w, int v);
    int getBit(int x, int w);
    double log2ceil(double x);

    int roulett();
    int findRandomDestination(int local_id,int hops);
    unsigned int getQueueSize() const;

    // Constructor
    SC_CTOR(ProcessingElement) {
	SC_METHOD(rxProcess);
	sensitive << reset;
	sensitive << clock.pos();

	SC_METHOD(txProcess);
	sensitive << reset;
	sensitive << clock.pos();
    }

};

#endif
