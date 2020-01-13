#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdexcept>
#include <sys/time.h>
#include <sstream>
#include <vector>

using namespace std;

const int max_chromosome_length = 1 << 28; //maximum length of a chromosome
const int max_lines_in_chromosome = 1 << 20; //maximum number of lines in a chromosome
const char encoding_rule[4] = {'A', 'C', 'G', 'T'};

int reference_sequence_len = 0;
int line_ending_array_len = 0;
int lowercase_array_len = 0;
int other_char_array_len = 0;
int letter_N_array_len = 0;

char *reference_sequence = new char[max_chromosome_length];
char *decompressed_sequence = new char[max_chromosome_length];
int *line_ending_indexes = new int[max_lines_in_chromosome];
int *line_ending_counts = new int[max_lines_in_chromosome];
int *lowercase_indexes = new int[max_lines_in_chromosome];
int *lowercase_counts = new int[max_lines_in_chromosome];
int *other_char_indexes = new int[max_lines_in_chromosome];
char *other_char_values = new char[max_lines_in_chromosome];
int *letter_N_indexes = new int[max_lines_in_chromosome];
int *letter_N_counts = new int[max_lines_in_chromosome];

vector<string> split(string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find(delimiter, pos_start)) != string::npos) {
        token = s.substr(pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back(token);
    }

    res.push_back(s.substr(pos_start));
    return res;
}

vector<string> split(char *chars, string delimiter) {
    string s(chars);
    return split(s, delimiter);
}

ifstream open_file_stream(char *file_path) {
    string ref_file(file_path);
    ifstream in("../" + ref_file);

    if (!in) {
        cout << "Unable to open file: " << ref_file << "\n";
        throw runtime_error("Cannot open file, closing app.");
    }
    return in;
}

/* Reads reference file and store's nucleotide bases' characters reference_sequence array */
void referenceFileToSequence(char *reference_file_path) {

    ifstream in = open_file_stream(reference_file_path);
    int str_len;
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
            if (temp_c == 'A' || temp_c == 'G' || temp_c == 'C' || temp_c == 'T') {
                reference_sequence[reference_sequence_len++] = temp_c;
            }
        }
    }
    reference_sequence[reference_sequence_len] = '\0';
    in.close();
}

void extractIndexAndCountInfo(ifstream *in, char *line, int *indexes, int *counts, int *counter) {
    (*in).getline(line, 2 << 19);
    vector<string> indexes_and_counts = split(line, " ");
    for (auto index_and_count_info : indexes_and_counts) {
        if (index_and_count_info == "") {
            break;
        }
        vector<string> index_and_count = split(index_and_count_info, "-");
        indexes[*counter] = stoi(index_and_count[0]);
        counts[(*counter)++] = stoi(index_and_count[1]);
    }
    for (int i = 0; i < *counter; i++) {
        cout << indexes[i] << '-' << counts[i] << '\n';
    }
}

void extractIndexAndCharInfo(ifstream *in, char *line, int *indexes, char *chars, int *counter) {
    (*in).getline(line, 2 << 19);
    vector<string> indexes_and_chars = split(line, " ");
    for (auto index_and_char_info : indexes_and_chars) {
        if (index_and_char_info == "") {
            break;
        }
        vector<string> index_and_char = split(index_and_char_info, "-");
        indexes[*counter] = stoi(index_and_char[0]);
        chars[(*counter)++] = index_and_char[1][0];
    }
    for (int i = 0; i < *counter; i++) {
        cout << indexes[i] << '-' << chars[i] << '\n';
    }
}

/* Extracts auxiliary data from target file */
void extractAuxiliaryInfoFromCompressedFile(char *filepath) {
    ifstream in = open_file_stream(filepath);
    char line[2 << 19];

    extractIndexAndCountInfo(&in, line, line_ending_indexes, line_ending_counts, &line_ending_array_len);
    extractIndexAndCountInfo(&in, line, lowercase_indexes, lowercase_counts, &lowercase_array_len);
    extractIndexAndCharInfo(&in, line, other_char_indexes, other_char_values, &other_char_array_len);
    extractIndexAndCountInfo(&in, line, letter_N_indexes, letter_N_counts, &letter_N_array_len);
}

void decompressToOutputFile(char *compressed_file_path, ofstream &output_file) {

    for (int i = 0; i < max_chromosome_length; i++) {
        decompressed_sequence[i] = '_';
    }

    for (int i = 0; i < other_char_array_len; i++) {
        decompressed_sequence[other_char_indexes[i]] = other_char_values[i];
    }

    for (int i = 0; i < letter_N_array_len; i++) {
        for (int j = letter_N_indexes[i]; j < letter_N_counts[i] + letter_N_indexes[i]; j++) {
            decompressed_sequence[j] = 'N';
        }
    }

    int current_index = 0;
    int j = 0;
    int output_sequence_len = 0;

    for (int i = 0; i < line_ending_array_len; i++) {
        int current_count = 0;
        while (current_count != line_ending_counts[i]) {
            current_index += line_ending_indexes[i];
            decompressed_sequence[current_index + j++] = '\n';
            current_count++;
        }
        output_sequence_len = current_index + j;
    }

    ifstream compressed_file = open_file_stream(compressed_file_path);
    char compressed_file_line[1024];

    int output_file_index = 0;
    int reference_sequence_index = 0;

    for (int i = 0; i < 5; i++) {
        compressed_file.getline(compressed_file_line, 1024);
    }

    while (compressed_file.getline(compressed_file_line, 1024)) {

        bool stop_condition = false;
        string mismatch_encoded_chars;
        bool matching = true;
        int index_of_match;
        int match_count;

        vector<string> index_and_count = split(compressed_file_line, " ");
        // if there is only one string after split, then we know these are mismatched characters
        if (index_and_count.size() == 1) {
            mismatch_encoded_chars = index_and_count[0];
            matching = false;
        } else {
            index_of_match = stoi(index_and_count[0]);
            match_count = stoi(index_and_count[1]);
            matching = true;
        }

        int length = matching ? match_count : mismatch_encoded_chars.size();

        for (int i = 0; i < length; i++) {
            if (decompressed_sequence[output_file_index] != '_') {
                output_file_index++;
                i--;
                continue;
            }
            char current_char;

            if (matching) {
                current_char = reference_sequence[index_of_match++];
            } else {
                current_char = encoding_rule[mismatch_encoded_chars[i] - '0'];
            }

            decompressed_sequence[output_file_index++] = current_char;
        }
    }

    for (int i = 0; i < lowercase_array_len; i++) {
        int index = lowercase_indexes[i];
        int count = lowercase_counts[i];
        for (int k = index; k < index + count; k++) {
            decompressed_sequence[k] = tolower(decompressed_sequence[k]);
        }
    }

    for (int i = 0; i < output_sequence_len; i++) {
        output_file << decompressed_sequence[i];
    }
}

int main(int argc, char *argv[]) {
    // init compressed and reference file pointers, output stream
    if (argc != 3) {
        cout << "Wrong number of arguments";
        return 0;
    }
    char *compressed_file = argv[1];
    char *reference_file = argv[2];
    ofstream output_file;
    output_file.open("decomp.txt");

    // init decompression timers
    struct timeval start;
    struct timeval end;
    unsigned long timer;
    gettimeofday(&start, NULL);

    //TODO insert 7za extraction

    // steps of the algorithm
    referenceFileToSequence(reference_file);
    extractAuxiliaryInfoFromCompressedFile(compressed_file);
    decompressToOutputFile(compressed_file, output_file);

    output_file.close();

    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("total decompression timer = %lf ms; %lf min\n", timer / 1000.0, timer / 1000.0 / 1000.0 / 60.0);

    return 0;
}
