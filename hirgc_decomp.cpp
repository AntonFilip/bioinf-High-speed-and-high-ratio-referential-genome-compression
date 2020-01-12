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
int reference_sequence_len = 0;
int line_ending_array_len = 0;
int lowercase_array_len = 0;
int other_char_len = 0;
int letter_N_array_len = 0;

char *reference_sequence = new char[max_chromosome_length];
char *target_sequence = new char[max_chromosome_length];
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

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
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
            if (temp_c == 'A' || temp_c == 'G' || temp_c == 'C' || temp_c == 'T' ) {
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
        cout << indexes[i] << '-' << counts[i] <<  '\n';
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
        cout << indexes[i] << '-' << chars[i] <<  '\n';
    }
}

/* Extracts auxiliary data from target file */
void extractAuxiliaryInfoFromCompressedFile(char *filepath) {
    ifstream in = open_file_stream(filepath);
    char line[2 << 19];

    extractIndexAndCountInfo(&in, line, line_ending_indexes, line_ending_counts, &line_ending_array_len);
    extractIndexAndCountInfo(&in, line, lowercase_indexes, lowercase_counts, &lowercase_array_len);
    extractIndexAndCharInfo(&in, line, other_char_indexes, other_char_values, &other_char_len);
    extractIndexAndCountInfo(&in, line, letter_N_indexes, letter_N_counts, &letter_N_array_len);
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

    output_file.close();

    gettimeofday(&end, NULL);
    timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
    printf("total decompression timer = %lf ms; %lf min\n", timer / 1000.0, timer / 1000.0 / 1000.0 / 60.0);

    return 0;
}
