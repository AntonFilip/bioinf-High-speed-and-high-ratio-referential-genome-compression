#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

using namespace std;

const int MAX_CHAR_NUM = 1 << 28; //maximum length of a chromosome
int *encoded_reference_sequence;
int *latest_index;
int *previous_index;
const int max_arr_num = 1<<30; // maximum length of hash table
int ref_sq_len;

int encoding_rule(char ch) {
    switch (ch) {
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

inline void init() {
    encoded_reference_sequence = new int[MAX_CHAR_NUM];
    for (int i = 0; i < max_arr_num; i++) {//initial entries
        latest_index[i] = -1;
    }
}

ifstream open_file_stream(char *arg) {
    string ref_file(arg);
    ifstream in("../" + ref_file);

    if (!in) {
        cout << "Unable to open file: " << ref_file << "\n";
        throw std::runtime_error("Cannot open file, closing app.");
    }
    return in;
}

/* Reads reference file and store's nucleotide bases' characters in 0-3 encoding to encoded_reference_sequence array */
void reference_file_to_encoded_sequence(char *arg) {

    ifstream in = open_file_stream(arg);
    int str_len, encoded;
    char temp_c;
    char str[1024];

    while (in.getline(str, 1024)) {
        str_len = strlen(str);
        for (int i = 0; i < str_len; i++) {
            temp_c = str[i];
            if (islower(temp_c)) {
                temp_c = toupper(temp_c);
            }
            encoded = encoding_rule(temp_c);
            if (encoded > -1) {
                encoded_reference_sequence[ref_sq_len++] = encoded;
            }
        }
    }
    in.close();
}

void read_target_file(char *arg) {

}

void construct_hash_table(int *encoded_reference_sequence) {

}

int main(int argc, char *argv[]) {

    init();

    //init reference and target file pointers
    char *ref_file = nullptr;
    char *target_file = nullptr;

    if (argc > 3 || argc < 2) {
        cout << "Wrong number of arguments";
        return 0;
    } else if (argc == 3) {
        target_file = argv[2];
    }
    ref_file = argv[1];

    reference_file_to_encoded_sequence(ref_file);
    construct_hash_table(encoded_reference_sequence);

    read_target_file(target_file);
    for (int i = 0; i < 1024; i++) {
        cout << encoded_reference_sequence[i] << " ";
    }
    return 0;
}
