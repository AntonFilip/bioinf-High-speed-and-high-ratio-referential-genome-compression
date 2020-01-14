# bioinf-High-speed-and-high-ratio-referential-genome-compression
Project for course Bioinformatics at University of Zagreb (Croatia), Faculty of Electrical Engineering and Computing https://www.fer.unizg.hr/predmet/bio referencing article: https://academic.oup.com/bioinformatics/article/33/21/3364/3885699.

Installation steps:

Clone repo
Compile compression source file with command: g++ -o hirgc hirgc_comp.cpp -std=c++11
Compile decompression source file with command: g++ -o decomp hirgc_decomp.cpp -std=c++11
To test compiled programs, you can use attached examples with commands:

./hirgc ce6_chrI.fa ce10_chrI.fa
./decomp ce10_chrI.fa_ref_ce6_chrI.fa ce6_chrI.fa
To ensure that both target and decompressed files are equal, use diff tool.
7za and diff tools should have executable permission that can be added using command:
chmod +x diff