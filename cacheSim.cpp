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


enum cache_line_status {MATCH, VACANCY, FULL}; //cache lookup results

typedef unsigned long int uli; //shortcut

class way{
	
	public:
	
	uli dirty;
	int tag;
	uli this_mod_time;

	way(){
		dirty=0;
		tag=-1;
		this_mod_time=0;
	}


};

class set_class{
	
	public:

	vector<way*> way_vec;
	uli last_mod_time;

	set_class(uli ways){
		for (int i; i<ways; i++){
			way_vec.push_back(new way);
		}

		last_mod_time=0;
	}

};


class cache_class{

	vector<set_class*> L1;
	vector<set_class*> L2;

	uli L1_cache_size_bytes;
	uli L2_cache_size_bytes;

	uli L1_set_num;
	uli L2_set_num;

	uli L1_ways;
	uli L2_ways;

	uli L1_index_bits;
	uli L2_index_bits;
		
	uli block_size_bytes;
	
	uli L1_access_time;
	uli L2_access_time;
	uli mem_access_time;

	uli L1_total_access;
	uli L2_total_access;

	uli L1_total_misses;
	uli L2_total_misses;

	uli total_access_time;
	uli write_allock;

	public:
	cache_class(uli block_size, uli L1_cache_size_bytes , uli L2_cache_size_bytes, 
				uli L1_assoc, uli L2_assoc, uli write_allock, uli L1_access_time, uli L2_access_time, uli mem_access_time){
					
		this->block_size_bytes = block_size;
		
		this->L1_cache_size_bytes = L1_cache_size_bytes;
		this->L2_cache_size_bytes = L2_cache_size_bytes;
		
		this->L1_ways = pow(2,L1_assoc);
		this->L2_ways = pow(2,L2_assoc);

		this->L1_set_num = L1_cache_size_bytes/(L1_ways*block_size_bytes); // calculate number of sets in each cache.
		this->L2_set_num = L2_cache_size_bytes/(L2_ways*block_size_bytes); // number of sets will be the length of the L1, L2 vectors.

		this->write_allock = write_allock;

		this->L1_index_bits = log2 (L1_set_num);
		this->L2_index_bits = log2 (L2_set_num);

		this-> L1_access_time    = L1_access_time;
		this-> L2_access_time    = L2_access_time;
		this-> mem_access_time   = mem_access_time;
		this-> total_access_time = 0;

		this->L1_total_access=0;
		this->L2_total_access=0;

		//init vectors
		for(int i=0; i<L1_set_num; i++){
			L1.push_back(new set_class(L1_ways));
		}

		for(int i=0; i<L2_set_num; i++){
			L2.push_back(new set_class(L2_ways));	
		}


	}

	uli get_L1_index(uli address){//find appropriate set in L1

		uli L1_index = address << 2;
		L1_index = L1_index<< int(log2(this->block_size_bytes));
		L1_index = L1_index % (L1_set_num);
		return L1_index;

	}

	uli get_L2_index(uli address){

		uli L2_index = address <<2;
		L2_index = L2_index << int(log2(this->block_size_bytes));
		L2_index = L2_index % (L2_set_num);
		return L2_index;

	}

	int get_LRU_way (set_class* set){
		int min_acc_index=0;

		for (int i =0; i<set->way_vec.size(); i++){
			if (set->way_vec[min_acc_index]->this_mod_time > set->way_vec[i]->this_mod_time){
				min_acc_index = i;
			}
		}
		return min_acc_index;
	}

	cache_line_status access_L1 (uli address, int* L1_way){

		uli L1_index = get_L1_index(address);
		int empty_block_flag =0;

		for(int i=0; i<L1_ways; i++){
			if(this->L1[L1_index]->way_vec[i]->tag==address){ //check if we already have this data in L1
				*L1_way = i;
				return MATCH;
			}
			if(this->L1[L1_index]->way_vec[i]->tag==-1){ //check if there is an empty way in this set
				*L1_way = i;
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

	cache_line_status access_L2 (uli address, int* L2_way, char command){
		uli L2_index = get_L2_index(address);
		int empty_block_flag =0;

		for(int i=0; i<L2_ways; i++){
			if(this->L2[L2_index]->way_vec[i]->tag==address){ //check if we already have this data in L1
				*L2_way = i;
				return MATCH;
			}
			if(this->L2[L2_index]->way_vec[i]->tag==-1){ //check if there is an empty way in this set
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


	void execute_command(char command, uli address){
		int L1_way =0, L2_way;
		uli L1_index = get_L1_index(address);
		cache_line_status L1_lookup_result = access_L1(address, &L1_way);
		this->L1_total_access ++;
		this->total_access_time += this->L1_access_time;

		if(L1_lookup_result == MATCH)//data found in L1
		{
			if (command=='w'){
				this->L1[L1_index]->way_vec[L1_way]->dirty = 1;
			}
			this->L1[L1_index]->last_mod_time += 1; //increase set max acces time by 1
			this->L1[L1_index]->way_vec[L1_way]->this_mod_time = this->L1[L1_index]->last_mod_time; //update block access time.
		}

		else{ //data not found in L1, trying L2
			this->L1_total_access ++;
			this->total_access_time += this->L2_access_time;
			cache_line_status L2_lookup_result = access_L2(address, &L2_way, command);
		}
		
		
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
