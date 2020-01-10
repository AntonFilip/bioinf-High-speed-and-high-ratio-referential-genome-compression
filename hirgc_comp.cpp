#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

using namespace std;

const int max_chromosome_length = 1 << 24; //maximum length of a chromosome
const int hash_table_length = 1 << 24; // maximum length of hash table
int encoded_reference_sequence_len = 0;
int k = 4;
int *encoded_reference_sequence = new int[max_chromosome_length];;
int *previous_index = new int[max_chromosome_length];
int *latest_index = new int[hash_table_length];

int *newline_indices = new int[1 << 20];
int *lower_sq_begin_indices = new int[1 << 20];
int *lower_sq_lengths = new int[1 << 20];
int *n_sq_begin_indices = new int[1 << 20];
int *n_sq_lengths = new int[1 << 20];
int *other_char_indices = new int[1 << 20];
char *other_chars = new char[1 << 20];

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
    for (int i = 0; i < hash_table_length; i++) {//initial entries
        latest_index[i] = -1;
    }
}

ifstream open_file_stream(char *file_path) {
    string ref_file(file_path);
    ifstream in("../" + ref_file);

    if (!in) {
        cout << "Unable to open file: " << ref_file << "\n";
        throw std::runtime_error("Cannot open file, closing app.");
    }
    return in;
}

/* Reads reference file and store's nucleotide bases' characters in 0-3 encoding to encoded_reference_sequence array */
void reference_file_to_encoded_sequence(char *reference_file_path) {

    ifstream in = open_file_stream(reference_file_path);
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
                encoded_reference_sequence[encoded_reference_sequence_len++] = encoded;
            }
        }
    }
    in.close();
}

/*creates auxiliary data from target file*/
void extract_auxiliary_info_from_tar_file(char *filepath) {
    ifstream in = open_file_stream(filepath);
    char str[1024];
    char id[100];
    char temp_char;
    int line_length;
    bool previous_upper = true;
    bool previous_n = true;
    int lower_sq_len = 0, n_sq_len = 0, other_char_len = 0, newline_len = 0;
    int file_pos_index = 0;

    in.getline(id, 100);

    while (in.getline(str, 1024)) {
        line_length = strlen(str);
        for (int i = 0; i < line_length; i++) {
            temp_char = str[i];
            if (islower(temp_char)) {
                if (previous_upper) {
                    previous_upper = false;
                    lower_sq_begin_indices[lower_sq_len] = file_pos_index;
                }
                temp_char = toupper(temp_char);
            }
            file_pos_index++;
        }
    }
}

/* constructs hash table from reference chromosome sequence */
void construct_hash_table(int *encoded_reference_sequence) {
    uint64_t tuple_value = 0; // number of bits in tuple_value has to be 2 * k
    for (int encoded_char_index = 0; encoded_char_index < encoded_reference_sequence_len; encoded_char_index++) {
        tuple_value = tuple_value << 2; // shift to left to make "room" for new encoded character
        tuple_value += encoded_reference_sequence[encoded_char_index]; // add new encoded character to tuple
        if (encoded_char_index < k - 1) { // used to skip first k - 1 values (characters) because we need exactly k values for tuples
            continue;
        }
        tuple_value = tuple_value & ((1 << 2 * k) - 1); // used to remove bits on indexes higher than 2 * k (older character that we don't need)
        int tuple_hash = tuple_value % hash_table_length; // tuple's hash is tuple's value modulo hash_table_length
        int tuple_index = encoded_char_index - (k - 1); //index of tuple if different from index of current character (different by k - 1)
        previous_index[tuple_index] = latest_index[tuple_hash];
        latest_index[tuple_hash] = tuple_index;
    }
}

void match_target_sequence_with_reference_and_output_to_file(FILE *resulting_file) {
    uint64_t tuple_value = 0; // number of bits in tuple_value has to be 2 * k
    for (int encoded_char_index = 0; encoded_char_index < encoded_reference_sequence_len; encoded_char_index++) {
        tuple_value = tuple_value << 2; // shift to left to make "room" for new encoded character
        tuple_value += encoded_reference_sequence[encoded_char_index]; // add new encoded character to tuple
        if (encoded_char_index < k - 1) { // used to skip first k - 1 values (characters) because we need exactly k values for tuples
            continue;
        }
        tuple_value = tuple_value & ((1 << 2 * k) - 1); // used to remove bits on indexes higher than 2 * k (older character that we don't need)
        int tuple_hash = tuple_value % hash_table_length; // tuple's hash is tuple's value modulo hash_table_length
        int tuple_index = encoded_char_index - (k - 1); //index of tuple if different from index of current character (different by k - 1)
        
}

int main(int argc, char *argv[]) {

    init();

    //init reference and target file pointers
    char *reference_file = nullptr;
    char *target_file = nullptr;

    if (argc > 3 || argc < 2) {
        cout << "Wrong number of arguments";
        return 0;
    } else if (argc == 3) {
        target_file = argv[2];
    }
    reference_file = argv[1];
    char resulting_file_name[100];
    sprintf(resulting_file_name, "%s_ref_%s", target_file, reference_file);
    FILE *resulting_file = fopen(resulting_file_name, "w");

    reference_file_to_encoded_sequence(reference_file);
    construct_hash_table(encoded_reference_sequence);
    extract_auxiliary_info_from_tar_file(target_file);
    match_target_sequence_with_reference_and_output_to_file(resulting_file);

    for (int i = 0; i < 1024; i++) {
        cout << encoded_reference_sequence[i] << " ";
    }
    return 0;
}
