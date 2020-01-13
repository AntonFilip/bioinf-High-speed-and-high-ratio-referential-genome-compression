#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
# include <sys/time.h>

using namespace std;

const int max_chromosome_length = 1 << 28; //maximum length of a chromosome
const int hash_table_length = 1 << 20; // maximum length of hash table
int encoded_reference_sequence_len = 0;
int encoded_target_sequence_len = 0;
int k = 20;
int newline_len = 0, lower_sq_len = 0, n_sq_len = 0, other_char_len = 0;

char id[100];
int *encoded_reference_sequence = new int[max_chromosome_length];;
int *previous_hashed_tuple_index = new int[max_chromosome_length];
int *last_hashed_tuple_index = new int[hash_table_length];

int *newline_indices = new int[max_chromosome_length];
int *lower_sq_begin_indices = new int[max_chromosome_length];
int *lower_sq_lengths = new int[max_chromosome_length];
int *n_sq_begin_indices = new int[max_chromosome_length];
int *n_sq_lengths = new int[max_chromosome_length];
int *other_char_indices = new int[max_chromosome_length];
char *other_chars = new char[max_chromosome_length];
int *encoded_target_sequence = new int[max_chromosome_length];

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
void referenceFileToEncodedSequence(char *reference_file_path) {

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

/* Extracts auxiliary data from target file */
void extractAuxiliaryAndSequenceInfoFromTargetFile(char *filepath) {
    ifstream in = open_file_stream(filepath);
    char str[1024];
    char temp_char;
    int line_length;
    bool previous_upper = true;
    bool previous_n = false;
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
        if (previous_n) {
            previous_n = false;
            n_sq_lengths[n_sq_len] = file_pos_index - n_sq_begin_indices[n_sq_len];
            n_sq_len++;
        }
        newline_indices[newline_len] = line_length;
        newline_len++;
        file_pos_index++;
    }

    if (!previous_upper) {
        lower_sq_lengths[lower_sq_len] = file_pos_index - lower_sq_begin_indices[lower_sq_len];
        lower_sq_len++;
    }

    if (previous_n) {
        n_sq_lengths[n_sq_len] = file_pos_index - n_sq_begin_indices[n_sq_len];
        n_sq_len++;
    }
}

/* Encodes data that has many repeating values by storing their value and running length.
 * Writes the encoded and formatted data to the given file.
 * Example: 21 21 21 21 21 21 5 5 5 7 7 7 7 7 7 will be encoded as 21-6 5-3 7-6*/

void encodeLineLengths(ofstream &outfile) {

    int first = newline_indices[0];
    int curr = newline_indices[0];
    int same_cnt = 0;

    for (int i = 1; i < newline_len + 1; i++) {

        if (curr == first) {
            same_cnt++;
        } else {
            outfile << first << "-" << same_cnt << " ";
            first = curr;
            same_cnt = 0;
            i--;
            continue;
        }
        curr = newline_indices[i];

    }
    outfile << first << "-" << same_cnt << "\n";

}

void outputTargetAuxiliaryData(ofstream &outfile) {

    //write lower sequence begin positions and length of those sequences
    int i = 0;
    for (i = 0; i < lower_sq_len - 1; i++) {
        outfile << lower_sq_begin_indices[i] << "-" << lower_sq_lengths[i] << " ";
    }
    outfile << lower_sq_begin_indices[i] << "-" << lower_sq_lengths[i] << "\n";

    //write other char positions and their values
    for (i = 0; i < other_char_len - 1; i++) {
        outfile << other_char_indices[i] << "-" << other_chars[i] << " ";
    }
    if (i > 0 || other_char_len == 1) {
        outfile << other_char_indices[i] << "-" << other_chars[i];
    }
    outfile << "\n";

    //write N sequence positions and their lengths
    for (i = 0; i < n_sq_len - 1; i++) {
        outfile << n_sq_begin_indices[i] << "-" << n_sq_lengths[i] << " ";
    }
    if (i > 0 || n_sq_len == 1) {
        outfile << n_sq_begin_indices[i] << "-" << n_sq_lengths[i];
    }
    outfile << "\n\n";
}

/* Constructs hash table from reference chromosome sequence */
void constructHashTable(int *encoded_reference_sequence) {
    // initialize last_hashed_tuple_index array
    for (int i = 0; i < hash_table_length; i++) {
        last_hashed_tuple_index[i] = -1;
    }

    uint64_t tuple_value = 0; // number of bits in tuple_value has to be 2 * k
    for (int encoded_char_index = 0; encoded_char_index < encoded_reference_sequence_len; encoded_char_index++) {
        tuple_value = tuple_value << 2; // shift to left to make "room" for new encoded character
        tuple_value += encoded_reference_sequence[encoded_char_index]; // add new encoded character to tuple
        if (encoded_char_index < k - 1) { // used to skip first k - 1 values (characters) because we need exactly k values for tuples
            continue;
        }
        if (k < 32) {
            tuple_value = tuple_value & (((uint64_t) 1 << (2 * k)) - 1); // used to remove bits on indexes higher than 2 * k (older character that we don't need)
        }
        int tuple_hash = tuple_value % (uint64_t) hash_table_length; // tuple's hash is tuple's value modulo hash_table_length
        int tuple_index = encoded_char_index - (k - 1); // index of tuple if different from index of current character (different by k - 1)
        previous_hashed_tuple_index[tuple_index] = last_hashed_tuple_index[tuple_hash];
        last_hashed_tuple_index[tuple_hash] = tuple_index;
    }
}

void matchTargetSequenceWithReferenceAndOutputToFile(ofstream &outfile) {
    int target_tuple_index = 0, mismatch_index = 0;
    while (target_tuple_index < encoded_target_sequence_len - k + 1) {

        uint64_t target_tuple_value = 0; // number of bits in target_tuple_value has to be 2 * k
        for (int j = 0; j < k; j++) {
            target_tuple_value = target_tuple_value << 2; // shift to left to make "room" for new encoded character
            target_tuple_value += encoded_target_sequence[target_tuple_index + j]; // add new encoded character to tuple
        }

        int target_tuple_hash = target_tuple_value % hash_table_length; // tuple's hash is tuple's value modulo hash_table_length
        int reference_tuple_index = last_hashed_tuple_index[target_tuple_hash]; // reference_tuple_index stores index of a tuple in reference sequence that has the same hash as target tuple

        int longest_match_index = -1;
        int longest_match_length = 0;

        while (reference_tuple_index != -1) { // repeat while there are more tuples left in the same bucket of a hash map
            int current_match_length = 0;
            // count how many characters match from target and reference sequences
            while ((reference_tuple_index + current_match_length < encoded_target_sequence_len) &&
                   (reference_tuple_index + current_match_length < encoded_reference_sequence_len) &&
                   encoded_reference_sequence[reference_tuple_index + current_match_length] ==
                   encoded_target_sequence[target_tuple_index + current_match_length]) {
                current_match_length += 1;
            }
            // if we have new longest match, store it
            if (current_match_length >= k && current_match_length > longest_match_length) {
                longest_match_index = reference_tuple_index;
                longest_match_length = current_match_length;
            }
            // fetch the next reference tuple (it's index) from a hash map bucket
            reference_tuple_index = previous_hashed_tuple_index[reference_tuple_index];
        }
        // if we found that target sequence matches some sequence in reference file, write that info to file and also write if there was any mismatch in between current and previous sequence
        if (longest_match_length > 0) {
            bool found_mismatch = false;
            for (int i = mismatch_index; i <= target_tuple_index - 1; i++) {
                found_mismatch = true;
                outfile << encoded_target_sequence[i];
            }
            if (found_mismatch) {
                outfile << "\n";
            }
            outfile << longest_match_index << " " << longest_match_length << "\n";
            mismatch_index = target_tuple_index + longest_match_length;
        }
        // jump to next sequence by adding longest match length to current index of target sequence
        target_tuple_index += longest_match_length + 1;
    }
    // write to file if there are any mismatches left at the end of the target file
    for (; mismatch_index < encoded_target_sequence_len; mismatch_index++) {
        outfile << encoded_target_sequence[mismatch_index];
    }
}

int main(int argc, char *argv[]) {
    // init reference and target file pointers, output stream
    if (argc != 3) {
        cout << "Wrong number of arguments";
        return 0;
    }

    char *reference_file = argv[1];
    char *target_file = argv[2];
    ofstream output_file;
    char terminal_command[150], output_file_name[100];

    sprintf(output_file_name, "%s_ref_%s", target_file, reference_file);

    output_file.open(output_file_name);

    // init compression timers
    struct timeval start;
    struct timeval end;
    unsigned long timer;
    gettimeofday(&start, NULL);

    // steps of the algorithm
    referenceFileToEncodedSequence(reference_file);
    constructHashTable(encoded_reference_sequence);
    extractAuxiliaryAndSequenceInfoFromTargetFile(target_file);
    //write chromosome identifier
    output_file << id << '\n';
    encodeLineLengths(output_file);
    outputTargetAuxiliaryData(output_file);
    matchTargetSequenceWithReferenceAndOutputToFile(output_file);

    output_file.close();

    sprintf(terminal_command, "./7za a %s.7z %s -m0=PPMd", output_file_name, output_file_name);
    system(terminal_command);
    sprintf(terminal_command, "rm %s", output_file_name);
    system(terminal_command);


    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("total compression timer = %lf ms; %lf min\n", timer / 1000.0, timer / 1000.0 / 1000.0 / 60.0);

    return 0;
}
