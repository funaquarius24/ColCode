/*
 * Noxim - the NoC Simulator 
 *
 * (C) 2005-2018 by the University of Catania
 * For the complete list of authors refer to file ../doc/AUTHORS.txt
 * For the license applied to these sources refer to file ../doc/LICENSE.txt
 *
 * This file contains the implementation of the processing element
 */

#include "ProcessingElement.h"

int ProcessingElement::randInt(int min, int max) {
    return min +
        (int)((double)(max - min + 1) * rand() / (RAND_MAX + 1.0));
}

void ProcessingElement::rxProcess() { 
    if (reset.read()) {
        ack_rx.write(0);
        current_level_rx = 0;
    } else {
        if (req_rx.read() == 1 - current_level_rx) {
            Flit flit_tmp = flit_rx.read();
            current_level_rx = 1 - current_level_rx; // Negate the old value for Alternating Bit Protocol (ABP)
            if(flit_tmp.flit_type == FLIT_TYPE_HEAD){
                clearDeps(flit_tmp);
            }
        }
        ack_rx.write(current_level_rx);
    }
}

void ProcessingElement::txProcess() {
    if (reset.read()) {
        req_tx.write(0);
        current_level_tx = 0;
        transmittedAtPreviousCycle = false;
    } else {
        Packet packet;

        if (canShot(packet)) {
            if(GlobalParams::traffic_distribution != TRAFFIC_MCSL)
                packet_queue.push(packet);
            transmittedAtPreviousCycle = true;
        } else
            transmittedAtPreviousCycle = false;

        if (ack_tx.read() == current_level_tx) {
            if (!packet_queue.empty()) {
                Flit flit = nextFlit(); // Generate a new flit
                flit_tx -> write(flit); // Send the generated flit
                current_level_tx = 1 - current_level_tx; // Negate the old value for Alternating Bit Protocol (ABP)
                req_tx.write(current_level_tx);
            }
        }
    }
}

Flit ProcessingElement::nextFlit() {
    Flit flit;
    Packet packet = packet_queue.front();
    //cout << "packet_size: " << packet.size << endl;

    flit.src_id = packet.src_id;
    flit.dst_id = packet.dst_id;
    flit.vc_id = packet.vc_id;
    flit.timestamp = packet.timestamp;
    flit.sequence_no = packet.size - packet.flit_left;
    flit.sequence_length = packet.size;
    flit.hop_no = 0;
    flit.task_edge = packet.task_edge;
    flit.packet_sequence = packet.packet_sequence;
    flit.iter_number = packet.iter_number;
    //  flit.payload     = DEFAULT_PAYLOAD;

    flit.hub_relay_node = NOT_VALID;

    if (packet.size == packet.flit_left)
        flit.flit_type = FLIT_TYPE_HEAD;
    else if (packet.flit_left == 1)
        flit.flit_type = FLIT_TYPE_TAIL;
    else
        flit.flit_type = FLIT_TYPE_BODY;

    packet_queue.front().flit_left--;
    if (packet_queue.front().flit_left == 0)
        packet_queue.pop();

    return flit;
}

bool ProcessingElement::canShot(Packet & packet) {
    //cout << "Before: " << endl;
    //assert(false);
    if (never_transmit) return false;

    //if(local_id!=16) return false;
    /* DEADLOCK TEST 
	double current_time = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;

	if (current_time >= 4100) 
	{
	    //if (current_time==3500)
	         //cout << name() << " IN CODA " << packet_queue.size() << endl;
	    return false;
	}
	//*/

    #ifdef DEADLOCK_AVOIDANCE
    if (local_id % 2 == 0)
        return false;
    #endif
    bool shot = false;
    double threshold;

    double now = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;

    if (GlobalParams::traffic_distribution != TRAFFIC_TABLE_BASED) {
        if(GlobalParams::traffic_distribution == TRAFFIC_MCSL )
        {
            if(now - GlobalParams::reset_time > GlobalParams::stats_warm_up_time){
                
                mcsl_ready = true;
                
                if(metDeps()){
                    if (!transmittedAtPreviousCycle)
                        threshold = GlobalParams::packet_injection_rate;
                    else
                        threshold = GlobalParams::probability_of_retransmission;
                    shot = (((double) rand()) / RAND_MAX < threshold);

                    if(shot){
                        if(list_of_taskInt.front().outgoing_edge_dep_list.size() > 0){
                            //cout << outgoing_dep_edge_list -> size() << endl;  
                            //packet = traffic_mcsl();
                            
                            traffic_mcsl();
                            
                            //assert(false);
                            if(!packets_left){
                                list_of_taskInt.front().outgoing_edge_dep_list.erase(list_of_taskInt.front().outgoing_edge_dep_list.begin());
                            }
                            
                            //outgoing_dep_edge_list -> erase(outgoing_dep_edge_list -> begin());
                            //cout << "after: " << endl;
                            return true;
                        }

                        // ism This means that there will be an extra cycle of no action. is this right?
                        list_of_taskInt.erase(list_of_taskInt.begin());
                        ongoing_task = false;
                    }    
                        
                }
            }
            
        }
        else{
            if (!transmittedAtPreviousCycle)
                threshold = GlobalParams::packet_injection_rate;
            else
                threshold = GlobalParams::probability_of_retransmission;

            shot = (((double) rand()) / RAND_MAX < threshold);
            if (shot) {
                if (GlobalParams::traffic_distribution == TRAFFIC_RANDOM)
                    packet = trafficRandom();
                else if (GlobalParams::traffic_distribution == TRAFFIC_TRANSPOSE1)
                    packet = trafficTranspose1();
                else if (GlobalParams::traffic_distribution == TRAFFIC_TRANSPOSE2)
                    packet = trafficTranspose2();
                else if (GlobalParams::traffic_distribution == TRAFFIC_BIT_REVERSAL)
                    packet = trafficBitReversal();
                else if (GlobalParams::traffic_distribution == TRAFFIC_SHUFFLE)
                    packet = trafficShuffle();
                else if (GlobalParams::traffic_distribution == TRAFFIC_BUTTERFLY)
                    packet = trafficButterfly();
                else if (GlobalParams::traffic_distribution == TRAFFIC_LOCAL)
                    packet = trafficLocal();
                else if (GlobalParams::traffic_distribution == TRAFFIC_ULOCAL)
                    packet = trafficULocal();
                else {
                    cout << "Invalid traffic distribution: " << GlobalParams::traffic_distribution << endl;
                    exit(-1);
                }
            }
        }
        
    } 
    else { // Table based communication traffic
        if (never_transmit)
            return false;

        bool use_pir = (transmittedAtPreviousCycle == false);
        vector < pair < int, double > > dst_prob;
        double threshold =
            traffic_table -> getCumulativePirPor(local_id, (int) now, use_pir, dst_prob);

        double prob = (double) rand() / RAND_MAX;
        shot = (prob < threshold);
        if (shot) {
            for (unsigned int i = 0; i < dst_prob.size(); i++) {
                if (prob < dst_prob[i].second) {
                    int vc = randInt(0, GlobalParams::n_virtual_channels - 1);
                    packet.make(local_id, dst_prob[i].first, vc, now, getRandomSize());
                    break;
                }
            }
        }
    }

    return shot;
}

void ProcessingElement::clearReadymadeDeps(){
    if(local_id == 53){
        //cout << "size: " << map_all_edge_packets_itr[2].size() << endl;
        //cout << "end" << endl;
    }//cout << "before: 224" << endl;

    if(!map_all_edge_packets_itr.count(mcslCount)) return;

    if(map_all_edge_packets_itr.at(mcslCount).size() > 0){
        //cout << "after: 246" << endl;
        for(uint i = 0; i < list_of_taskInt.size(); i++){
            int j = 0;
            vector<pair<int, int>>::iterator it;
            //cout << "before: ";
            //cout << "after: 234 " << endl;
            //cout << list_of_taskInt.at(i).incoming_edge_dep_list.size() << endl;
            //int te = 0;
            for(it = list_of_taskInt.at(i).incoming_edge_dep_list.begin(); 
                it != list_of_taskInt.at(i).incoming_edge_dep_list.end(); it++, j++){
                //cout << te++ << endl;
                int edgeID = it->first; //cout << it << endl;
                map<int, pair<int, int>>::iterator edge = map_all_edge_packets_itr[mcslCount].find(edgeID);
                if(edge != map_all_edge_packets_itr[mcslCount].end()){
                    if(edge->second.first == edge->second.second){
                        list_of_taskInt.at(i).incoming_edge_dep_list.erase(it);
                        it--;
                    }
                    
                }
            }//cout << "after: " << endl;
            //cout << "after: 247" << endl;
        }//cout << "after: 248" << endl;
        
    }
        

}


bool ProcessingElement::metDeps(){
    //cout << local_id << endl;
    //cout << "Before: ";
    if(!ongoing_task){
        //cout << "after: 256" << endl;
        startTaskTime = sc_time_stamp().to_double();
        if(finishedTaskIteration){
            //cout << "REACHED!!!" << endl;
            
            mcslCount++;
            int numb = 300;
            if(mcslCount >= numb) return false;
            if(mcslCount >= numb -1){
                //cout << "ID: " << local_id << " mcslCount: " << mcslCount << endl;
                //cout << local_id << ", " ;
            }
            local_rec_proc = *rec_proc;

            

            taskList = &local_rec_proc.task_list;
            //local_rec_proc.task_list_pos = 0;
            
            for(uint i = 0; i < taskList -> size(); i++) {
                TaskInt t;
                t.id = taskList -> at(i) -> get_id();
                t.position = i;
                rec_task = taskList -> at(i);
                for(uint j = 0; j < rec_task->incoming_edge_list.size(); j++){
                    t.incoming_edge_dep_list.push_back(pair<int, int>(rec_task->incoming_edge_list.at(j)->get_id(), j));
                }
                for(uint j = 0; j < rec_task->outgoing_edge_list.size(); j++){
                    t.outgoing_edge_dep_list.push_back(pair<int, int>(rec_task->outgoing_edge_list.at(j)->get_id(), j));
                }
                list_of_taskInt.push_back(t);
                
                //task_list_int.push_back(pair<int, int>(taskList -> at(i) -> get_id(), i));
            }//cout << "after: 288" << endl;
            clearReadymadeDeps();
           // << "after: 292" << endl;
            
            //cout << "localID: " << local_id << endl;
            
            
        }

        if(list_of_taskInt.size() == 0){
            finishedTaskIteration = true;
            if(mcslCount == 1){
                //cout << local_id << ", " ;
            }

            map_all_edge_packets_itr.erase((mcslCount) ); 

            return false;
        }
        
        //rec_task = taskList -> at(task_list_int.front().second);
        rec_task = taskList -> at(list_of_taskInt.front().position);

        execTimes = &rec_task -> get_recorded_execution_time();
        execTime = execTimes -> at(mcslCount % 20); 

        if(local_id > -1){
            //cout << "mcslCount: " << mcslCount << " incoming: " << list_of_taskInt.front().incoming_edge_dep_list.size() << " ID: " << local_id <<  endl;
        }

        ongoing_task = true;
        finishedTaskIteration = false;
        
        //assert(false);
    }
    else{
        int now =  sc_time_stamp().to_double();
        if(now - startTaskTime >= execTime && list_of_taskInt.front().incoming_edge_dep_list.size() == 0){
            
            return true;
        }//cout << "after: "<< endl;
    }
    //cout << "after: here" << endl;
    return false;
    
}

void ProcessingElement::clearDeps(Flit &flit){
    //cout << "Before: " << endl;
    //if(local_id == 3)
        //cout << "PB: " << local_id << " mcslCount: " << mcslCount << " task: " << rec_task->get_id() << endl;
    int flit_task_id = flit.task_edge.first;
    //RecTask* task = local_rec_proc.get_task(flit_task_id);

    if(local_id == 2){
        //cout << "ID: " << local_id << "mcslCount: " <<  flit.iter_number << " " << flit.task_edge.first << " count: " << count++ << endl;
    }

    //map_all_edge_packets[flit.task_edge.second]++;
    map_all_edge_packets_itr[flit.iter_number][flit.task_edge.second].first++; 
    map_all_edge_packets_itr[flit.iter_number][flit.task_edge.second].second = flit.packet_sequence.first; 
    if(flit.iter_number != mcslCount){
        return;
    }

    //if(map_all_edge_packets[flit.task_edge.second] == flit.packet_sequence.first){
    if(map_all_edge_packets_itr[flit.iter_number][flit.task_edge.second].first == flit.packet_sequence.first){ 
        
        for(uint i = 0; i < list_of_taskInt.size(); i++){
            TaskInt *t = &list_of_taskInt.at(i);
            if(flit.task_edge.first == t ->id) {
                for(uint j = 0; j < t->incoming_edge_dep_list.size(); j++){
                    if(flit.task_edge.second == t->incoming_edge_dep_list.at(j).first){
                        t->incoming_edge_dep_list.erase(t->incoming_edge_dep_list.begin() + j);
                        //cout << "after: " << endl;
                        return;
                    }
                }

            }
            if(local_id == 2){
                //cout << mcslCount << " taskId: " << flit_task_id << " numRecv: " << map_all_edge_packets_itr[flit.iter_number][flit.task_edge.second];
                //cout << " packLength: " << flit.packet_sequence.first << " task List size: " << task_list_int.size() << endl;
            }  
            //cout << "Before: " << endl;
            //cout << "mcslCount: " << mcslCount << ( mcslCount > flit.iter_number ) << endl;
            //if(flit.task_edge.second == task_list_int.at(i).first && flit.iter_number == mcslCount){
                //task -> incoming_edge_list.erase(task -> incoming_edge_list.begin() + i);
                //cout << "Before: " << endl;
                
                //cout << "After: " << endl;
            //}//cout << "After: " << endl;
        }
    }//cout << "after: here" << endl;
    
    
    
}

// ism This generates the MCSL traffic
void ProcessingElement::traffic_mcsl(){

    if(!packets_left){
        
        
        rec_edge = *rec_task->outgoing_edge_list.at(list_of_taskInt.front().outgoing_edge_dep_list.front().second);
        int dataSize = ceil( rec_edge.get_recorded_msg_size().at(mcslCount % 20) );
        //cout << "PB: " << local_id << " edge: " << rec_edge.get_id() << " task: " << rec_task -> get_id() << " dataSize: " << dataSize << endl;
        
        
        // ism Since i'm dealing with words, 1 word is 16 bits
        // how many words make one flit?
        int words = GlobalParams::flit_size / 16;
        int num_flits = dataSize / words;
        
        int num_packets = ceil((num_flits + 0.0) / GlobalParams::max_packet_size);
        double timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;

        

        
        //Packet p;
        packet_tmpl.src_id = local_id; 
        
        packet_tmpl.dst_id = rec_edge.get_dst_proc_id();
        packet_tmpl.task_edge = pair<int, int>(rec_edge.get_dst_task_id(), rec_edge.get_id());
        
        packet_tmpl.packet_sequence = pair<int, int>(num_packets, 0);

        //packet_tmpl.timestamp = timestamp;

        packet_tmpl.iter_number = mcslCount;

        
        packets_left = true;
            packet_number = 0;
        
    }else
    {
        packet_number++;
        Packet p;
        p.src_id = packet_tmpl.src_id;
        
        p.dst_id = packet_tmpl.dst_id;
        p.task_edge = pair<int, int>(packet_tmpl.task_edge.first, packet_tmpl.task_edge.second);
        
        p.packet_sequence = pair<int, int>(packet_tmpl.packet_sequence.first, packet_number);
        

        p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
        p.size = p.flit_left = GlobalParams::max_packet_size;
        p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);

        p.iter_number = packet_tmpl.iter_number;

        
        packet_queue.push(p);

        

        if(packet_number == p.packet_sequence.first){
            packets_left = false;
            //packet_number = 0;
        }
        if(local_id == 2){
            //cout << sc_time_stamp().to_double()/1000 << " Packet queued!" << endl;
        }
    }
    

    /* //cout << "before: " << endl;
    //rec_edge = *outgoing_dep_edge_list -> at(outgoing_dep_edge_list_int.front().second);
    //rec_edge = rec_task.outgoing_edge_list -> at(list_of_taskInt.front().outgoing_edge_dep_list.front().second);
    rec_edge = *rec_task->outgoing_edge_list.at(list_of_taskInt.front().outgoing_edge_dep_list.front().second);
    int dataSize = ceil( rec_edge.get_recorded_msg_size().at(mcslCount % 20) );
    //cout << "PB: " << local_id << " edge: " << rec_edge.get_id() << " task: " << rec_task -> get_id() << " dataSize: " << dataSize << endl;
    
    
    // ism Since i'm dealing with words, 1 word is 16 bits
    // how many words make one flit?
    int words = GlobalParams::flit_size / 16;
    int num_flits = dataSize / words;
    
    int num_packets = ceil((num_flits + 0.0) / GlobalParams::max_packet_size);
    double timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;

    

    for(int i = 0; i < num_packets; i++){
        Packet p;
        p.src_id = local_id;
        
        p.dst_id = rec_edge.get_dst_proc_id();
        p.task_edge = pair<int, int>(rec_edge.get_dst_task_id(), rec_edge.get_id());
        
        p.packet_sequence = pair<int, int>(num_packets, i + 1);

        p.timestamp = timestamp;
        p.size = p.flit_left = GlobalParams::max_packet_size;
        p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);

        p.iter_number = mcslCount;

        
        packet_queue.push(p);
        
    }
     */
    //cout << "after" << endl;
    return;

}

Packet ProcessingElement::trafficLocal() {
    Packet p;
    p.src_id = local_id;
    p.iter_number = mcslCount;
    double rnd = rand() / (double) RAND_MAX;

    vector < int > dst_set;

    int max_id = (GlobalParams::mesh_dim_x * GlobalParams::mesh_dim_y);

    for (int i = 0; i < max_id; i++) {
        if (rnd <= GlobalParams::locality) {
            if (local_id != i && sameRadioHub(local_id, i))
                dst_set.push_back(i);
        } else
        if (!sameRadioHub(local_id, i))
            dst_set.push_back(i);
    }

    int i_rnd = rand() % dst_set.size();

    p.dst_id = dst_set[i_rnd];
    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();
    p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);

    return p;
}

int ProcessingElement::findRandomDestination(int id, int hops) {
    assert(GlobalParams::topology == TOPOLOGY_MESH);

    int inc_y = rand() % 2 ? -1 : 1;
    int inc_x = rand() % 2 ? -1 : 1;

    Coord current = id2Coord(id);

    for (int h = 0; h < hops; h++) {

        if (current.x == 0)
            if (inc_x < 0) inc_x = 0;

        if (current.x == GlobalParams::mesh_dim_x - 1)
            if (inc_x > 0) inc_x = 0;

        if (current.y == 0)
            if (inc_y < 0) inc_y = 0;

        if (current.y == GlobalParams::mesh_dim_y - 1)
            if (inc_y > 0) inc_y = 0;

        if (rand() % 2)
            current.x += inc_x;
        else
            current.y += inc_y;
    }
    return coord2Id(current);
}

int roulette() {
    int slices = GlobalParams::mesh_dim_x + GlobalParams::mesh_dim_y - 2;

    double r = rand() / (double) RAND_MAX;

    for (int i = 1; i <= slices; i++) {
        if (r < (1 - 1 / double(2 << i))) {
            return i;
        }
    }
    assert(false);
    return 1;
}

Packet ProcessingElement::trafficULocal() {
    Packet p;
    p.src_id = local_id;

    int target_hops = roulette();

    p.dst_id = findRandomDestination(local_id, target_hops);

    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();

    return p;
}

Packet ProcessingElement::trafficRandom() {
    Packet p;
    p.src_id = local_id;
    double rnd = rand() / (double) RAND_MAX;
    double range_start = 0.0;
    int max_id;

    if (GlobalParams::topology == TOPOLOGY_MESH)
        max_id = (GlobalParams::mesh_dim_x * GlobalParams::mesh_dim_y) - 1; //Mesh 
    else // other delta topologies
        max_id = GlobalParams::n_delta_tiles - 1;

    // Random destination distribution
    do {
        p.dst_id = randInt(0, max_id);

        // check for hotspot destination
        for (size_t i = 0; i < GlobalParams::hotspots.size(); i++) {

            if (rnd >= range_start && rnd < range_start + GlobalParams::hotspots[i].second) {
                if (local_id != GlobalParams::hotspots[i].first) {
                    p.dst_id = GlobalParams::hotspots[i].first;
                }
                break;
            } else
                range_start += GlobalParams::hotspots[i].second; // try next
        }
        #ifdef DEADLOCK_AVOIDANCE
        assert((GlobalParams::topology == TOPOLOGY_MESH));
        if (p.dst_id % 2 != 0) {
            p.dst_id = (p.dst_id + 1) % 256;
        }
        #endif

    } while (p.dst_id == p.src_id);

    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();
    p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);

    return p;
}
// TODO: for testing only
Packet ProcessingElement::trafficTest() {
    Packet p;
    p.src_id = local_id;
    p.dst_id = 10;

    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();
    p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);

    return p;
}

Packet ProcessingElement::trafficTranspose1() {
    assert(GlobalParams::topology == TOPOLOGY_MESH);
    Packet p;
    p.src_id = local_id;
    Coord src, dst;

    // Transpose 1 destination distribution
    src.x = id2Coord(p.src_id).x;
    src.y = id2Coord(p.src_id).y;
    dst.x = GlobalParams::mesh_dim_x - 1 - src.y;
    dst.y = GlobalParams::mesh_dim_y - 1 - src.x;
    fixRanges(src, dst);
    p.dst_id = coord2Id(dst);

    p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);
    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();

    return p;
}

Packet ProcessingElement::trafficTranspose2() {
    assert(GlobalParams::topology == TOPOLOGY_MESH);
    Packet p;
    p.src_id = local_id;
    Coord src, dst;

    // Transpose 2 destination distribution
    src.x = id2Coord(p.src_id).x;
    src.y = id2Coord(p.src_id).y;
    dst.x = src.y;
    dst.y = src.x;
    fixRanges(src, dst);
    p.dst_id = coord2Id(dst);

    p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);
    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();

    return p;
}

void ProcessingElement::setBit(int & x, int w, int v) {
    int mask = 1 << w;

    if (v == 1)
        x = x | mask;
    else if (v == 0)
        x = x & ~mask;
    else
        assert(false);
}

int ProcessingElement::getBit(int x, int w) {
    return (x >> w) & 1;
}

inline double ProcessingElement::log2ceil(double x) {
    return ceil(log(x) / log(2.0));
}

Packet ProcessingElement::trafficBitReversal() {

    int nbits =
        (int)
    log2ceil((double)
        (GlobalParams::mesh_dim_x *
            GlobalParams::mesh_dim_y));
    int dnode = 0;
    for (int i = 0; i < nbits; i++)
        setBit(dnode, i, getBit(local_id, nbits - i - 1));

    Packet p;
    p.src_id = local_id;
    p.dst_id = dnode;

    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();

    return p;
}

Packet ProcessingElement::trafficShuffle() {

    int nbits =
        (int)
    log2ceil((double)
        (GlobalParams::mesh_dim_x *
            GlobalParams::mesh_dim_y));
    int dnode = 0;
    for (int i = 0; i < nbits - 1; i++)
        setBit(dnode, i + 1, getBit(local_id, i));
    setBit(dnode, 0, getBit(local_id, nbits - 1));

    Packet p;
    p.src_id = local_id;
    p.dst_id = dnode;

    p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);
    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();

    return p;
}

Packet ProcessingElement::trafficButterfly() {

    int nbits = (int) log2ceil((double)
        (GlobalParams::mesh_dim_x *
            GlobalParams::mesh_dim_y));
    int dnode = 0;
    for (int i = 1; i < nbits - 1; i++)
        setBit(dnode, i, getBit(local_id, i));
    setBit(dnode, 0, getBit(local_id, nbits - 1));
    setBit(dnode, nbits - 1, getBit(local_id, 0));

    Packet p;
    p.src_id = local_id;
    p.dst_id = dnode;

    p.vc_id = randInt(0, GlobalParams::n_virtual_channels - 1);
    p.timestamp = sc_time_stamp().to_double() / GlobalParams::clock_period_ps;
    p.size = p.flit_left = getRandomSize();

    return p;
}

void ProcessingElement::fixRanges(const Coord src,
    Coord & dst) {
    // Fix ranges
    if (dst.x < 0)
        dst.x = 0;
    if (dst.y < 0)
        dst.y = 0;
    if (dst.x >= GlobalParams::mesh_dim_x)
        dst.x = GlobalParams::mesh_dim_x - 1;
    if (dst.y >= GlobalParams::mesh_dim_y)
        dst.y = GlobalParams::mesh_dim_y - 1;
}

int ProcessingElement::getRandomSize() {
    return randInt(GlobalParams::min_packet_size,
        GlobalParams::max_packet_size);
}


unsigned int ProcessingElement::getQueueSize() const {
    return packet_queue.size();
}