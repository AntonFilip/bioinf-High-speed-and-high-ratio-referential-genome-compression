#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>

using namespace std;

const int MAX_CHAR_NUM = 1<<28;//maximum length of a chromosome
int *ref_sq_tuple_values;
int ref_sq_len;

int encoding_rule(char ch) {
    switch(ch){
        case 'A':
            return 0;
        case 'C':
            return 1;
        case 'G':
            return 2;
        case 'T':
            return 3;
        default:
            return -1;
    }
}

inline void init(){
    ref_sq_tuple_values = new int[MAX_CHAR_NUM];
}

void read_reference_file(char *arg){

    string ref_file(arg);
    ifstream in("../" + ref_file);

    if(!in) {
        cout << "Unable to open file: " << ref_file << "\n";
        return;
    }

    int str_len, index;
    char temp_c;
    char str[1024];

    while(in.getline(str,1024)) {
        str_len = strlen(str);
        for(int i = 0; i < str_len; i++){
            temp_c = str[i];
            if(islower(temp_c)){
                temp_c = toupper(temp_c);
            }
            index = encoding_rule(temp_c);
            if (index > -1){
                ref_sq_tuple_values[ref_sq_len++] = index;
            }
        }
    }
    in.close();
}

void read_target_file(char *arg){

}

int main(int argc, char *argv[]) {

    init();
    char *ref_file = nullptr;
    char *target_file = nullptr;

    if(argc > 3 || argc < 2){
        cout << "Wrong number of arguments";
        return 0;
    } else if(argc == 3) {
        target_file = argv[2];
    }
    ref_file = argv[1];
    read_reference_file(ref_file);
    for(int i = 0; i < 1024; i++){
        cout << ref_sq_tuple_values[i] << " ";
    }
    return 0;
}
