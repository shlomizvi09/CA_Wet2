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
	

	uli write_allock;

	public:
	cache_class(uli block_size, uli L1_cache_size_bytes , uli L2_cache_size_bytes, 
				uli L1_assoc, uli L2_assoc, uli write_allock){
					
		block_size_bytes = block_size;
		
		L1_cache_size_bytes = L1_cache_size_bytes;
		L2_cache_size_bytes = L2_cache_size_bytes;
		
		L1_ways = pow(2,L1_assoc);
		L2_ways = pow(2,L2_assoc);

		L1_set_num = L1_cache_size_bytes/(L1_ways*block_size_bytes); // calculate number of sets in each cache.
		L2_set_num = L2_cache_size_bytes/(L2_ways*block_size_bytes); // number of sets will be the length of the L1, L2 vectors.

		write_allock = write_allock;

		L1_index_bits = log2 (L1_set_num);
		L2_index_bits = log2 (L2_set_num);

		//init vectors
		for(int i=0; i<L1_set_num; i++){
			L1.push_back(new set_class(L1_ways));
		}

		for(int i=0; i<L2_set_num; i++){
			L2.push_back(new set_class(L2_ways));	
		}


	}

	uli get_L1_index(uli address){

		uli L1_index = address <<2;
		L1_index = address << int(log2(this->block_size_bytes));
		L1_index = L1_index % (L1_ways);
		return L1_index;

	}

	uli get_L2_index(uli address){

		uli L2_index = address <<2;
		L2_index = address << int(log2(this->block_size_bytes));
		L2_index = L2_index % (L2_ways);
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

		//TODO: add mem access time here

		if(empty_block_flag){
			
			this->L2[L2_index]->way_vec[*L2_way]->tag=address;
			return VACANCY;
		}
		else{
			
			if((this->write_allock && command=='w')||command=='r'){ //evict and replace LRU way

				*L2_way = get_LRU_way(this->L2[L2_index]);
				this->L2[L2_index]->way_vec[*L2_way]->tag=address;
				this->L2[L2_index]->last_mod_time +=1; //increase set max acces time by 1
				this->L2[L2_index]->way_vec[*L2_way]->this_mod_time = this->L2[L2_index]->last_mod_time; //update block access time.

				if( command=='w'){
					this->L2[L2_index]->way_vec[*L2_way]->dirty=1; //mark block as dirty
				}
				else{
					this->L2[L2_index]->way_vec[*L2_way]->dirty=0;
				}
				
			}
			return FULL;
		}

	 }


	void execute_command(char command, uli address){
		int L1_way =0, L2_way;
		uli L1_index = get_L1_index(address);
		cache_line_status L1_lookup_result = access_L1(address, &L1_way);

		switch (L1_lookup_result)
		{
		case MATCH:
			if (command=='w'){
				L1[L1_index]->way_vec[L1_way]->dirty = 1;
			}
			this->L1[L1_index]->last_mod_time +=1; //increase set max acces time by 1
			this->L1[L1_index]->way_vec[L1_way]->this_mod_time = this->L1[L1_index]->last_mod_time; //update block access time.
			//TODO: add L2 access time calculation here
		break;

		case VACANCY:
		//TODO: add  L2 access time calculation here
			access_L2(address, &L2_way, command);
		break;

		case FULL:

			
			if((this->write_allock && command=='w')||command=='r'){ //evict and replace LRU way in cache L1

				L1_way = get_LRU_way(this->L1[L1_index]);
				this->L1[L1_index]->way_vec[L1_way]->tag=address;
				this->L1[L1_index]->last_mod_time +=1; //increase set max acces time by 1
				this->L1[L1_index]->way_vec[L1_way]->this_mod_time = this->L1[L1_index]->last_mod_time; //update block access time.

				if( command=='w'){
					this->L1[L1_index]->way_vec[L1_way]->dirty=1; //mark block as dirty
				}
				else{
					this->L1[L1_index]->way_vec[L1_way]->dirty=0;
				}
				
			}
			break;

		default:
			assert(0 && "should not reach here");
			break;
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
	cache_class cache(BSize,L1Size, L2Size, L1Assoc, L2Assoc,WrAlloc);

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
