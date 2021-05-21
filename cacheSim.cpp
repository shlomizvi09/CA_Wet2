/* 046267 Computer Architecture - Spring 2021 - HW #2 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <cassert>

#define ADDRESS_LENGTH 30

using std::vector;
using std::FILE;
using std::string;
using std::cout;
using std::endl;
using std::cerr;
using std::ifstream;
using std::stringstream;


enum cache_line_status {DEFAULT=-1, MATCH, VACANCY, FULL}; //cache lookup results

typedef uint32_t uint; //shortcut

class way{
	
	public:
	
	uint dirty;
	int tag;
	uint this_mod_time;

	way(){
		dirty=0;
		tag=-1;
		this_mod_time=0;
	}

	bool is_dirty(){
		return dirty == 0 ? false : true;
	}

	void set_tag(uint tag){
		if (tag == 0){
			cerr << "tag is zero";
			return;
		}
		this->tag = tag;
	 	return;
	}

	uint get_tag(){
		return this->tag;
	}


};

class set_class{
	
	public:

	vector<way> way_vec;
	uint last_mod_time;

	set_class(uint ways) : way_vec(ways){
		assert(way_vec.size() == ways);
		last_mod_time=0;
	}

	int get_way_index(uint address){
		for (std::vector<way>::size_type i = 0; i < this->way_vec.size(); i++){
			if (this->way_vec[i].tag == address){
				return i;
			}
		}
		return -1;
	}

	int get_LRU_way_index (){
		int min_acc_index=0;

		for (std::vector<way>::size_type i = 0; i < this->way_vec.size(); i++){
			if (this->way_vec[min_acc_index].this_mod_time > this->way_vec[i].this_mod_time){
				min_acc_index = i;
			}
		}

		return min_acc_index;
	}

};


class cache_class{
	vector<set_class> L1;
	vector<set_class> L2;

	uint L1_cache_size_bytes;
	uint L2_cache_size_bytes;

	uint L1_set_num;
	uint L2_set_num;

	uint L1_ways;
	uint L2_ways;

	uint L1_index_bits;
	uint L2_index_bits;
		
	uint block_size_bytes;
	
	uint L1_access_time;
	uint L2_access_time;
	uint mem_access_time;

	uint L1_total_access;
	uint L2_total_access;

	uint L1_total_misses;
	uint L2_total_misses;

	uint total_access_time;
	uint write_alloc;

	public:
	cache_class(uint block_size, uint L1_cache_size_bytes , uint L2_cache_size_bytes, 
				uint L1_assoc, uint L2_assoc, uint write_alloc, uint L1_access_time, uint L2_access_time,
				uint mem_access_time) : L1(L1_set_num), L2(L2_set_num){

		this->block_size_bytes = block_size;
		
		this->L1_cache_size_bytes = L1_cache_size_bytes;
		this->L2_cache_size_bytes = L2_cache_size_bytes;
		
		this->L1_ways = pow(2,L1_assoc);
		this->L2_ways = pow(2,L2_assoc);

		this->L1_set_num = L1_cache_size_bytes/(L1_ways*block_size_bytes); // calculate number of sets in each cache.
		this->L2_set_num = L2_cache_size_bytes/(L2_ways*block_size_bytes); // number of sets will be the length of the L1, L2 vectors.

		this->write_alloc = write_alloc;

		this->L1_index_bits = log2 (L1_set_num);
		this->L2_index_bits = log2 (L2_set_num);

		this-> L1_access_time    = L1_access_time;
		this-> L2_access_time    = L2_access_time;
		this-> mem_access_time   = mem_access_time;
		this-> total_access_time = 0;

		this->L1_total_access=0;
		this->L2_total_access=0;
		this->L1_total_misses=0;
		this->L2_total_misses=0;

		this->write_alloc = write_alloc;

		this->L1_index_bits = log2 (L1_set_num);
		this->L2_index_bits = log2 (L2_set_num);


	}

	uint get_L1_set_index(uint address){//find appropriate set in L1

		uint L1_index = address << 2;
		L1_index = L1_index<< int(log2(this->block_size_bytes));
		L1_index = L1_index % (L1_set_num);
		return L1_index;

	}

	uint get_L2_set_index(uint address){

		uint L2_index = address <<2;
		L2_index = L2_index << int(log2(this->block_size_bytes));
		L2_index = L2_index % (L2_set_num);
		return L2_index;

	}

	double get_L1_miss_rate(){
		return this->L1_total_misses / this->L1_total_access;
	}

	double get_L2_miss_rate(){
		return this->L2_total_misses / this->L2_total_access;
	}

	cache_line_status access_L1 (uint address, uint* L1_way_index){

		uint L1_index = get_L1_set_index(address);
		int empty_block_flag =0;

		for(int i=0; i<L1_ways; i++){
			if(this->L1[L1_index].way_vec[i].tag==address){ //check if we already have this data in L1
				*L1_way_index = i;
				return MATCH;
			}
			if(this->L1[L1_index].way_vec[i].tag==-1){ //check if there is an empty way in this set
				*L1_way_index = i;
				empty_block_flag=1;
			}
		}
		if(empty_block_flag){
			return VACANCY;
		}
		else{
			return FULL;
		}
	 }

	cache_line_status access_L2 (uint address, uint* L2_way, char command){
		uint L2_index = get_L2_set_index(address);
		int empty_block_flag =0;

		for(int i=0; i<L2_ways; i++){
			if(this->L2[L2_index].way_vec[i].tag==address){ //check if we already have this data in L1
				*L2_way = i;
				return MATCH;
			}
			if(this->L2[L2_index].way_vec[i].tag==-1){ //check if there is an empty way in this set
				*L2_way = i;
				empty_block_flag=1;
			}
		}
		if(empty_block_flag){
			return VACANCY;
		}
		else{
			return FULL;
		}
	 }

	void evict_block_from_L1(uint L1_set_index, uint L1_way_index){
		assert(this->L1[L1_set_index].way_vec[L1_way_index].tag != -1);
		if (this->L1[L1_set_index].way_vec[L1_way_index].dirty){  // Write Back
			uint L2_set_index = this->get_L2_set_index(this->L1[L1_set_index].way_vec[L1_way_index].tag);
			uint L2_way_index = this->L2[L2_set_index].get_way_index(this->L1[L1_set_index].way_vec[L1_way_index].tag);
			this->L2[L2_set_index].way_vec[L2_way_index].dirty = 1;
		}
		this->L1[L1_set_index].way_vec[L1_way_index].dirty = 0;
		this->L1[L1_set_index].way_vec[L1_way_index].tag = -1;
	}

	void evict_block_from_L2(uint L2_set_index, uint L2_way_index){
		uint address = this->L2[L2_set_index].way_vec[L2_way_index].tag;
		assert(address != -1);
		uint L1_way_index = 0;
		cache_line_status L1_lookup_res = access_L1(address, &L1_way_index);
		if (L1_lookup_res == MATCH){
			uint L1_set_index = this->get_L1_set_index(address);
			this->evict_block_from_L1(L1_set_index, L1_way_index);
		}
		if (this->L2[L2_set_index].way_vec[L2_way_index].dirty){
			// do nothing, it's a Write Back.
		}
		this->L2[L2_set_index].way_vec[L2_way_index].dirty = 0;
		this->L2[L2_set_index].way_vec[L2_way_index].tag = -1;
	}

	// void insert_block_to_L1(uint address){

	// }

	void execute_command(char command, uint address){
		uint L1_way_index = 0, L2_way_index = 0;
		uint L1_set_index = get_L1_set_index(address);
		uint L2_set_index = this->get_L2_set_index(address);
		cache_line_status L1_lookup_result = access_L1(address, &L1_way_index);
		cache_line_status L2_lookup_result = DEFAULT;
		this->L1_total_access ++;
		this->total_access_time += this->L1_access_time;

		if(L1_lookup_result == MATCH)//data found in L1
		{
			if (command=='w'){
				this->L1[L1_set_index].way_vec[L1_way_index].dirty = 1;
			}
			this->L1[L1_set_index].last_mod_time += 1; //increase set max access time by 1
			this->L1[L1_set_index].way_vec[L1_way_index].this_mod_time = this->L1[L1_set_index].last_mod_time; //update block access time.
			return;
		}
		this->L1_total_misses++;
		L2_lookup_result = access_L2(address, &L2_way_index, command); //data not found in L1, trying L2
		this->L2_total_access ++;
		this->total_access_time += this->L2_access_time;
		if (L2_lookup_result == MATCH){
			if (command == 'w' && !(this->write_alloc)){ // Writing in "No Write Allocate" method, without copying to L1
				this->L2[L2_set_index].way_vec[L2_way_index].dirty = 1;
				this->L2[L2_set_index].last_mod_time += 1; //increase set max access time by 1
				this->L2[L2_set_index].way_vec[L2_way_index].this_mod_time = this->L2[L2_set_index].last_mod_time; //update block access time.
				return;
			}
			// if we got here, command is 'READ' or command is 'WRITE' with "Write Alloc", so we bring the block to L1 anyway.
			if (L1_lookup_result = FULL){
				uint L1_LRU_way_index = this->L1[L1_set_index].get_LRU_way_index();
				L1_way_index = L1_LRU_way_index;
				this->evict_block_from_L1(L1_set_index, L1_LRU_way_index);
			}
			this->L1[L1_set_index].way_vec[L1_way_index].tag = address;
			this->L1[L1_set_index].way_vec[L1_way_index].dirty = (command == 'w') ? 1 : 0; // the difference between 'READ' and "Write-Alloc".
			this->L1[L1_set_index].last_mod_time += 1; //increase set max access time by 1
			this->L1[L1_set_index].way_vec[L1_way_index].this_mod_time = this->L1[L1_set_index].last_mod_time; //update block access time.
			return;
		}
		this->L2_total_misses++;
		// if we got here, the block isn't in L2 as well, so we access Mem.
		this->total_access_time += this->mem_access_time;
		if (command == 'w' && !(this->write_alloc)){
			return; // we don't actually write in this simulator. suppose we wrote, and return.
		}
		if (L2_lookup_result == FULL){
			uint L2_LRU_way_index = this->L2[L2_set_index].get_LRU_way_index();
			L2_way_index = L2_LRU_way_index;
			this->evict_block_from_L2(L2_set_index, L2_LRU_way_index);
		}
		this->L2[L2_set_index].way_vec[L2_way_index].tag = address;
		this->L2[L2_set_index].way_vec[L2_way_index].dirty = 0; // the difference between 'READ' and "Write-Alloc".
		this->L2[L2_set_index].last_mod_time += 1; //increase set max access time by 1
		this->L2[L2_set_index].way_vec[L2_way_index].this_mod_time = this->L2[L2_set_index].last_mod_time; //update block access time.
		L1_lookup_result = access_L1(address, &L1_way_index); //TODO: think if need to add L1 access time
		if (L1_lookup_result = FULL){
			uint L1_LRU_way_index = this->L1[L1_set_index].get_LRU_way_index();
			L1_way_index = L1_LRU_way_index;
			this->evict_block_from_L1(L1_set_index, L1_LRU_way_index);
		}
		this->L1[L1_set_index].way_vec[L1_way_index].tag = address;
		this->L1[L1_set_index].way_vec[L1_way_index].dirty = (command == 'w') ? 1 : 0; // the difference between 'READ' and "Write-Alloc".
		this->L1[L1_set_index].last_mod_time += 1; //increase set max access time by 1
		this->L1[L1_set_index].way_vec[L1_way_index].this_mod_time = this->L1[L1_set_index].last_mod_time; //update block access time.
		return;
	}

};

int main(int argc, char **argv) {

	if (argc < 19) {
		cerr << "Not enough arguments" << endl;
		return 0;
	}

	// Get input arguments

	// File
	// Assuming it is the first argument
	char* fileString = argv[1];
	ifstream file(fileString); //input file stream
	string line;
	if (!file || !file.good()) {
		// File doesn't exist or some other error
		cerr << "File not found" << endl;
		return 0;
	}

	unsigned MemCyc = 0, BSize = 0, L1Size = 0, L2Size = 0, L1Assoc = 0,
			L2Assoc = 0, L1Cyc = 0, L2Cyc = 0, WrAlloc = 0;

	for (int i = 2; i < 19; i += 2) {
		string s(argv[i]);
		if (s == "--mem-cyc") {
			MemCyc = atoi(argv[i + 1]);
		} else if (s == "--bsize") {
			BSize = atoi(argv[i + 1]);
		} else if (s == "--l1-size") {
			L1Size = atoi(argv[i + 1]);
		} else if (s == "--l2-size") {
			L2Size = atoi(argv[i + 1]);
		} else if (s == "--l1-cyc") {
			L1Cyc = atoi(argv[i + 1]);
		} else if (s == "--l2-cyc") {
			L2Cyc = atoi(argv[i + 1]);
		} else if (s == "--l1-assoc") {
			L1Assoc = atoi(argv[i + 1]);
		} else if (s == "--l2-assoc") {
			L2Assoc = atoi(argv[i + 1]);
		} else if (s == "--wr-alloc") {
			WrAlloc = atoi(argv[i + 1]);
		} else {
			cerr << "Error in arguments" << endl;
			return 0;
		}
	}
	cache_class cache(BSize,L1Size, L2Size, L1Assoc, L2Assoc,WrAlloc, L1Cyc, L2Cyc, MemCyc);

	while (getline(file, line)) {

		stringstream ss(line);
		string address;
		char operation = 0; // read (R) or write (W)
		if (!(ss >> operation >> address)) {
			// Operation appears in an Invalid format
			cout << "Command Format error" << endl;
			return 0;
		}

		// DEBUG - remove this line
		cout << "operation: " << operation;

		string cutAddress = address.substr(2); // Removing the "0x" part of the address

		// DEBUG - remove this line
		cout << ", address (hex)" << cutAddress;

		unsigned long int num = 0;
		num = strtoul(cutAddress.c_str(), NULL, 16);

		// DEBUG - remove this line
		cout << " (dec) " << num << endl;

		cache.execute_command(operation, num);

	}

	double L1MissRate;
	double L2MissRate;
	double avgAccTime;

	printf("L1miss=%.03f ", L1MissRate);
	printf("L2miss=%.03f ", L2MissRate);
	printf("AccTimeAvg=%.03f\n", avgAccTime);

	return 0;
}
