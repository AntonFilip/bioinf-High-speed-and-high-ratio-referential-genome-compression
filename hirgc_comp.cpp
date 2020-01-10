#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>

using namespace std;

const int max_chromosome_length = 1 << 24; //maximum length of a chromosome
const int hash_table_length = 1 << 20; // maximum length of hash table
int encoded_reference_sequence_len = 0;
int encoded_target_sequence_len = 0;
int k = 20;
int *encoded_reference_sequence = new int[max_chromosome_length];;
int *previous_hashed_tuple_index = new int[max_chromosome_length];
int *last_hashed_tuple_index = new int[hash_table_length];

int *newline_indices = new int[1 << 20];
int *lower_sq_begin_indices = new int[1 << 20];
int *lower_sq_lengths = new int[1 << 20];
int *n_sq_begin_indices = new int[1 << 20];
int *n_sq_lengths = new int[1 << 20];
int *other_char_indices = new int[1 << 20];
char *other_chars = new char[1 << 20];
int *encoded_target_sequence = new int[1 << 20];

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
        last_hashed_tuple_index[i] = -1;
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
    in.getline(str, 1024);

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
    bool previous_n = false;
    int lower_sq_len = 0, n_sq_len = 0, other_char_len = 0, newline_len = 0;
    int file_pos_index = 0;
    int current_encoded_char = -1;

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
            } else {
                if (!previous_upper) {
                    previous_upper = true;
                    lower_sq_lengths[lower_sq_len] = file_pos_index - lower_sq_begin_indices[lower_sq_len];
                    lower_sq_len++;
                }
            }


            if (!previous_n) {
                if (temp_char == 'N') {
                    previous_n = true;
                    n_sq_begin_indices[n_sq_len] = file_pos_index;
                }
            } else {
                if (temp_char != 'N') {
                    previous_n = false;
                    n_sq_lengths[n_sq_len] = file_pos_index - n_sq_begin_indices[n_sq_len];
                    n_sq_len++;
                }
            }

            if (temp_char != 'N') {
                current_encoded_char = encoding_rule(temp_char);
                if (current_encoded_char > -1) {
                    encoded_target_sequence[encoded_target_sequence_len++] = current_encoded_char;
                } else {
                    other_char_indices[other_char_len] = file_pos_index;
                    other_chars[other_char_len] = temp_char;
                    other_char_len++;
                }
            }

            file_pos_index++;
        }
        newline_indices[newline_len] = line_length;
        newline_len++;
    }

    if (!previous_upper) {
        lower_sq_lengths[lower_sq_len] = file_pos_index - lower_sq_begin_indices[lower_sq_len];
        lower_sq_len++;
    }

    if (previous_n) {
        n_sq_lengths[n_sq_len] = file_pos_index - n_sq_begin_indices[n_sq_len];
        n_sq_len++;
    }

//    for (int i = 0; i < encoded_target_sequence_len; i++) {
//        cout << encoded_target_sequence[i];
//    }
//
//    cout << "\n";
//    cout << "\n";
//
//    for (int i = 0; i < lower_sq_len; i++) {
//        cout << lower_sq_begin_indices[i] << " " << lower_sq_lengths[i] << "\n";
//    }
//
//    cout << "\n";
//    cout << "\n";
//
//    for (int i = 0; i < other_char_len; i++) {
//        cout << other_char_indices[i] << " " << other_chars[i] << "\n";
//    }
//
//    cout << "\n";
//    cout << "\n";
//
//    for (int i = 0; i < n_sq_len; i++) {
//        cout << n_sq_begin_indices[i] << " " << n_sq_lengths[i] << "\n";
//    }
//
//    cout << "\n";
//    cout << "\n";
//
//    for (int i = 0; i < newline_len; i++) {
//        cout << newline_indices[i] << " " ;
//    }
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
        if (k < 32) {
            tuple_value = tuple_value & (((uint64_t)1 << (2 * k)) - 1); // used to remove bits on indexes higher than 2 * k (older character that we don't need)
        }
        int tuple_hash = tuple_value % hash_table_length; // tuple's hash is tuple's value modulo hash_table_length
        int tuple_index = encoded_char_index - (k - 1); //index of tuple if different from index of current character (different by k - 1)
        previous_hashed_tuple_index[tuple_index] = last_hashed_tuple_index[tuple_hash];
        last_hashed_tuple_index[tuple_hash] = tuple_index;
    }
}

void match_target_sequence_with_reference_and_output_to_file(FILE *resulting_file) {
    int target_tuple_index, mismatch_index = 0;
    while (target_tuple_index < encoded_target_sequence_len - k + 1) {
        uint64_t target_tuple_value = 0; // number of bits in target_tuple_value has to be 2 * k
        for (int j = 0; j < k; j++) {
            target_tuple_value = target_tuple_value << 2; // shift to left to make "room" for new encoded character
            target_tuple_value += encoded_target_sequence[target_tuple_index + j]; // add new encoded character to tuple
        }
        int target_tuple_hash = target_tuple_value % hash_table_length; // tuple's hash is tuple's value modulo hash_table_length

        int reference_tuple_index = last_hashed_tuple_index[target_tuple_hash];

        int longest_match_index = -1;
        int longest_match_length = 0;

        while (reference_tuple_index != -1) {
            int current_match_length = 0;
            while ((reference_tuple_index + current_match_length < encoded_target_sequence_len) && (reference_tuple_index + current_match_length < encoded_reference_sequence_len) && encoded_reference_sequence[reference_tuple_index + current_match_length] == encoded_target_sequence[target_tuple_index + current_match_length]){
                current_match_length += 1;
            }
            if (current_match_length > k && current_match_length > longest_match_length) {
                longest_match_index = reference_tuple_index;
                longest_match_length = current_match_length;
            }
            reference_tuple_index = previous_hashed_tuple_index[reference_tuple_index];
        }
        if (longest_match_length > 0) {
            bool foundMismatch = false;
            for (int i = mismatch_index; i <= target_tuple_index - 1; i++) {
                foundMismatch = true;
                fprintf(resulting_file, "%d", encoded_target_sequence[i]);
            }
            if (foundMismatch) {
                fprintf(resulting_file, "\n");
            }
            fprintf(resulting_file, "%d %d\n", longest_match_index, longest_match_length);
            mismatch_index = target_tuple_index + longest_match_length;
        }
        target_tuple_index += longest_match_length + 1;
    }
    for(; target_tuple_index < encoded_target_sequence_len; target_tuple_index++) {
        fprintf(resulting_file, "%d", encoded_target_sequence[target_tuple_index]);
    }
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
    FILE *resulting_file = fopen("result.txt", "w");

    reference_file_to_encoded_sequence(reference_file);
    construct_hash_table(encoded_reference_sequence);
    extract_auxiliary_info_from_tar_file(target_file);
    match_target_sequence_with_reference_and_output_to_file(resulting_file);

    fclose(resulting_file);
    return 0;
}
