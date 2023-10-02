Included HCTree.cpp, HCTree.hpp, Helper.cpp, Helper.hpp, MakeFile, Compress.cpp, Uncompress.cpp.

## HCTree.cpp
```
#include <queue>
#include <vector>
#include <fstream>
#include "Helper.hpp"
#include "HCTree.hpp"
#include <algorithm>
using namespace std;


void helpdelete(HCNode* root){
    if(root == NULL){
        return;
    }
    //recursion
    if(root->c0 != NULL){
        helpdelete(root->c0);
    }
    //recursion
    if(root->c1 != NULL){
        helpdelete(root->c1);
    }
    delete root;
}

HCTree::~HCTree(){
    helpdelete(root);
}


void HCTree::build(const vector<int>& freqs){
    //initialize
    priority_queue<HCNode*, vector<HCNode*>, HCNodePtrComp> pq;
    //construct HCNode
    //PRECONDITION: freqs is a vector of ints, such that freqs[i] is the frequency of occurrence of byte i in the input file.
    for(int i = 0; i < (int)freqs.size(); i++){
        int count = freqs[i];
        unsigned char symbol = i;
    
        //if count = 0, means won't happen
        if(count != 0){
            HCNode* newnode = new HCNode(count, symbol);
            pq.push(newnode);
            //and leaves[i] points to the leaf node containing byte i.
            leaves[i] = newnode;
        }
    }
    //loop if not only one tree
    while(pq.size() != 1){
        //pop fist one
        HCNode* left = pq.top(); pq.pop();
        //pop second one
        HCNode* right = pq.top();pq.pop();
        //combine 
        int sum = left->count + right->count;
        unsigned char newsymbol = right->symbol; // get the bigger one symbol
        //construct
        HCNode* cmb = new HCNode(sum, newsymbol);
        cmb->c0 = left;
        cmb->c1 = right;
        left->p = cmb;
        right->p = cmb;
        pq.push(cmb);
    }
    //POSTCONDITION: root points to the root of the trie,
    root = pq.top();
}


void HCTree::encode(unsigned char symbol, FancyOutputStream & out) const{
    HCNode* codingnode = leaves[symbol];
    //initialize
    vector<char> coding;
    HCNode* curr = codingnode;
    //more then one
    while(curr != root){
        HCNode* currpar = curr->p;
        //left child
        if(currpar->c0 == curr){
            coding.push_back('0');
        }
        //right child
        else if(currpar->c1 == curr){
            coding.push_back('1');
        }
        curr = currpar;
    }
    //empty tree
    if(root == nullptr){
        return;
    }
    //tree only root
    if(root->c0 == nullptr && root->c1 == nullptr){
        out.write_bit(0);
        return;
    }
    //output
    for(int i = coding.size()-1; i >= 0 ; i--){
        if(coding[i] == '0'){
            out.write_bit(0);
        }
        else if(coding[i] == '1'){
            out.write_bit(1);
        }
    }
}


unsigned char HCTree::decode(FancyInputStream & in) const{
    HCNode* curr = root;
    //empty
    if(root == nullptr){
        return 0;
    }
    //tree only root
    if(root->c0 == nullptr && root->c1 == nullptr){
        in.read_bit();
        return root->symbol;
    }
    //loop until break
    while(1){
        int code = in.read_bit();
        if(code == 0){
            curr = curr->c0;
        }
        else if(code == 1){
            curr = curr->c1;
        }
        //break: left node found
        if(curr->c0 == nullptr && curr->c1 == nullptr){
            break;
        }
    }
    return curr->symbol;
} 

HCNode* HCTree::helpbuild(FancyInputStream& in, HCNode* root){
    while(1){
        int bit = in.read_bit();
        if(root->c0 != nullptr && root->c1 != nullptr){
            return nullptr;
        }
        if(bit == 0){
            if(root != nullptr){
                root = new HCNode(0,0);
                helpbuild(in, root);
            }else{
                HCNode* tmp = new HCNode(0,0);
                tmp->p = root;
                tmp = root->c0 == nullptr ? root->c0 : root->c1;
                helpbuild(in, tmp);
                }
        }else{
            unsigned char bt = in.read_byte();
            HCNode* leaf = new HCNode(0, bt);
            if(root == nullptr){
                return leaf;
            }
            leaves[bt] = leaf;
            leaf->p = root;
            leaf = root->c0 == 0? root->c0 : root->c1;
        }
    }
}

int HCTree::helpdec(FancyInputStream& in){
    int count = in.read_int();
    helpbuild(in, root);
    return count;
}
```

## HCTree.hpp
```
#ifndef HCTREE_HPP
#define HCTREE_HPP
#include <queue>
#include <vector>
#include <fstream>
#include "Helper.hpp"
using namespace std;

/**
 * A Huffman Code Tree class
 */
class HCTree {
    public:
        HCNode* root;
        vector<HCNode*> leaves;

    public:
        /**
         * Constructor, which initializes everything to null pointers
         */
        HCTree() : root(nullptr) {
            leaves = vector<HCNode*>(256, nullptr);
        }

        ~HCTree();

        /**
         * Use the Huffman algorithm to build a Huffman coding tree.
         * PRECONDITION: freqs is a vector of ints, such that freqs[i] is the frequency of occurrence of byte i in the input file.
         * POSTCONDITION: root points to the root of the trie, and leaves[i] points to the leaf node containing byte i.
         */
        void build(const vector<int>& freqs);

        /**
         * Write to the given FancyOutputStream the sequence of bits coding the given symbol.
         * PRECONDITION: build() has been called, to create the coding tree, and initialize root pointer and leaves vector.
         */
        void encode(unsigned char symbol, FancyOutputStream & out) const;

        /**
         * Return symbol coded in the next sequence of bits from the stream.
         * PRECONDITION: build() has been called, to create the coding tree, and initialize root pointer and leaves vector.
         */
        unsigned char decode(FancyInputStream & in) const;

        HCNode* helpbuild(FancyInputStream& in, HCNode* root);

        int helpdec(FancyInputStream& in);
};
#endif // HCTREE_HPP
```

## Helper.cpp

```
#include "Helper.hpp"

// error function implementation
void error(const char* message) {
    cerr << "ERROR: " << message << endl; exit(1);
}

// FancyInputStream function implementations
FancyInputStream::FancyInputStream(const char* filename) : FILENAME(filename), input_file(ifstream(filename,ios::binary)), buffer(0), buffer_index(8) {}

bool FancyInputStream::good() const {
    return input_file.good();
}

int FancyInputStream::filesize() const  {
    return ifstream(FILENAME, ios::ate | ios::binary).tellg();
}

void FancyInputStream::reset() {
    input_file.clear();  // clear EOF flag
    input_file.seekg(0); // move to begining of file
    buffer = 0;          // clear bitwise buffer
    buffer_index = 8;    // move bitwise buffer index back to beginning
}

int FancyInputStream::read_int() {
    if(buffer_index != 8) {
        error("Attempting to read int when bitwise buffer is not empty");
    }
    unsigned int num;                          // temporary variable to store the number
    input_file.read((char*)&num, sizeof(num)); // read the number and store it in 'num'
    if(input_file.eof()) {                     // not enough bytes in the file to read it
        error("Not enough bytes to read the next int");
    }
    return num;
}

int FancyInputStream::read_byte() {
    return input_file.get();
}

int FancyInputStream::read_bit() {
    // if there are no more bits to read in the buffer,
    if(buffer_index == 8) {
        int const & tmp = read_byte(); // try to read the next byte

        // if there are no more bytes to read, there are no more bits to read
        if(tmp == -1) {
            return -1;
        }

        // we read a byte successfully, so update our buffer
        buffer = tmp;
        buffer_index = 0;
    }

    // read the next bit from the bitwise buffer
    return (buffer >> (7-buffer_index++)) & 1;
}

// FancyOutputStream function implementations
FancyOutputStream::FancyOutputStream(const char* filename) : output_file(ofstream(filename,ios::binary)), buffer(0), buffer_index(0) {}

FancyOutputStream::~FancyOutputStream() {
    flush();
}

bool FancyOutputStream::good() const {
    return output_file.good();
}

void FancyOutputStream::write_int(int const & num) {
    if(buffer_index != 0) {
        error("Attempting to write int when bitwise buffer is not empty");
    }
    output_file.write((char*)&num, sizeof(num)); // write 'num' to file
}

void FancyOutputStream::write_byte(unsigned char const & byte) {
    if(buffer_index != 0) {
        error("Attempting to write byte when bitwise buffer is not empty");
    }
    output_file.put(byte);
}

void FancyOutputStream::write_bit(int bit) {
    // crash if invalid input
    if(bit != 0 && bit != 1) {
        error("Trying to write invalid bit");
    }

    // add bit to bitwise buffer
    buffer |= (bit << (7-buffer_index++));

    // if the bitwise buffer is full,
    if(buffer_index == 8) {
        flush_bitwise(); // flush it
    }
}

void FancyOutputStream::flush_bitwise() {
    // if we have bits in our bitwise buffer,
    if(buffer_index != 0) {
        output_file.put(buffer); // write the bitwise buffer to the ofstream
        buffer = 0;              // reset the buffer
        buffer_index = 0;        // reset the buffer index
    }
}

void FancyOutputStream::flush() {
    flush_bitwise();     // try to flush the bitwise buffer
    output_file.flush(); // flush the ofstream
}

// HCNode function implementations
HCNode::HCNode(int count, unsigned char symbol) : count(count), symbol(symbol), c0(nullptr), c1(nullptr), p(nullptr) {}

bool HCNode::operator<(const HCNode& other) {
    // if the counts are different, compare counts
    if(count != other.count){
        return count > other.count;
    }

    // if the counts are equal, use symbol to break tie
    return symbol > other.symbol;
}

bool HCNodePtrComp::operator()(HCNode*& lhs, HCNode*& rhs) const {
    return *lhs < *rhs;
}
```

## Helper.hpp
```
#ifndef HELPER
#define HELPER
#include <fstream>
#include <iostream>
using namespace std;

/**
 * Conveniently crash with error messages
 */
void error(const char* message);

/**
 * Handle reading from a file. You must never call read_byte or read_int after calling read_bit!!!
 */
class FancyInputStream {
    private:
        // member variables (aka instance variables)
        const char* FILENAME; // input file's name
        ifstream input_file;  // input stream from which to read
        unsigned char buffer; // bitwise buffer (holds 8 bits = 1 byte)
        int buffer_index;     // current index of bitwise buffer

    public:
        /**
         * Constructor, which initializes a FancyInputStream object to read from the given file
         */
        FancyInputStream(const char* filename);

        /**
         * Returns true if none of the stream's error state flags (eofbit, failbit and badbit) is set.
         * See: https://www.cplusplus.com/reference/ios/ios/good/
         */
        bool good() const;

        /**
         * Return the size of the input file
         */
        int filesize() const;

        /**
         * Move back to the beginning of the input file and clear bitwise buffer
         */
        void reset();

        /**
         * Read a single (usually 4-byte) integer from the file,
         * or crash if there are not enough bytes left in the file
         */
        int read_int();

        /**
         * Read a single byte from the file as an int in the range [0,255],
         * or return -1 if we've reached the end of the file and no more bytes exist to read
         */
        int read_byte();

        /**
         * Read a single bit from the file as an int that is either 0 or 1,
         * or return -1 if we've reached the end of the file and no more bits exist to read
         */
        int read_bit();
};

/**
 * Handle writing to a file. You must never call write_byte or write_int after calling write_bit!!!
 */
class FancyOutputStream {
    public:
        // member variables (aka instance variables)
        ofstream output_file; // output stream to which to write
        unsigned char buffer; // bitwise buffer (holds 8 bits = 1 byte)
        int buffer_index;     // current index of bitwise buffer

    public:
        /**
         * Constructor, which initializes a FancyOutputStream object to write to the given file
         */
        FancyOutputStream(const char* filename);

        /**
         * Destructor, which flushes everything
         */
        ~FancyOutputStream();

        /**
         * Returns true if none of the stream's error state flags (eofbit, failbit and badbit) is set.
         * See: https://www.cplusplus.com/reference/ios/ios/good/
         */
        bool good() const;

        /**
         * Write a single (usually 4-byte) integer to the file
         */
        void write_int(int const & num);

        /**
         * Write a single byte to the file
         */
        void write_byte(unsigned char const & byte);

        /**
         * Write a single bit to the file
         */
        void write_bit(int bit);

        /**
         * Flush the bitwise buffer to the ofstream
         */
        void flush_bitwise();

        /**
         * Flush everything to the file
         */
        void flush();
};

/**
 * Represent nodes in an HCTree (Huffman Tree) object
 */
class HCNode {
    public:
        // member variables (aka instance variables)
        int count;            // count of this node
        unsigned char symbol; // symbol of this node
        HCNode* c0;           // pointer to '0' child
        HCNode* c1;           // pointer to '1' child
        HCNode* p;            // pointer to parent

        /**
         * Constructor, which initializes an HCNode object with a given count and symbol
         */
        HCNode(int count, unsigned char symbol);

        /**
         * Less-than operator to compare HCNodes deterministically
         */
        bool operator<(const HCNode& other);
};

/**
 * A 'function class' for use as the Compare class in a priority_queue<HCNode*>.
 * For this to work, operator< must be defined to do the right thing on HCNodes.
 */
class HCNodePtrComp {
    public:
        bool operator()(HCNode*& lhs, HCNode*& rhs) const;
};
#endif // HELPER
```

## Make File
```
# use g++ with C++11 support
CXX=g++
CXXFLAGS?=-Wall -pedantic -g -O0 -std=c++11
OUTFILES=compress uncompress

all: $(OUTFILES)

compress: compress.cpp Helper.cpp HCTree.cpp Helper.hpp HCTree.hpp
	$(CXX) $(CXXFLAGS) -o compress compress.cpp Helper.cpp HCTree.cpp

uncompress: uncompress.cpp Helper.cpp HCTree.cpp Helper.hpp HCTree.hpp
	$(CXX) $(CXXFLAGS) -o uncompress uncompress.cpp Helper.cpp HCTree.cpp

clean:
	rm -f $(OUTFILES) *.o

```

## Compress.cpp
```
#include <queue>
#include <vector>
#include <fstream>
#include <iostream>
#include "Helper.hpp"
#include "HCTree.hpp"
#include <algorithm>
using namespace std;


void helpencode(FancyOutputStream& out, HCNode* root, int count){
    if(root == nullptr){
        return;
    }
    if(count != -1){
        out.write_int(count);
    }
    if(root->c0 != nullptr || root->c1 != nullptr){
        out.write_bit('0');
    }else{
        out.write_bit('1');
        out.write_byte(root->symbol);
    }
    out.write_int('0');
    helpencode(out, root->c0, -1);
    helpencode(out, root->c1, -1);
}

int main(int argc, char* argv[]){
    
    //check number of parameters
    if(argc != 3){
        error("invalid number of arguments");
    }
    ifstream infile(argv[1]);
    ofstream outfile(argv[2]);
    //input file empty
    FancyInputStream intmp(argv[1]);
    if(intmp.filesize() == 0){
        return 0;
    }
    vector<int> freqs(256, 0);
    int count = 0;
    unsigned char next;
    next = intmp.read_byte();
    while(intmp.good() == true){
        freqs[next]++;
        next = intmp.read_byte();
        count++;
    }
    //build
    HCTree* tree = new HCTree();
    tree->build(freqs);

    FancyOutputStream fcout(argv[2]);
    FancyOutputStream& output_file = fcout;
    helpencode(output_file, tree->root, count);
    while(intmp.good()){
        next = intmp.read_byte();
        tree->encode(next, output_file);
    }
    if(output_file.buffer_index != 0){
        output_file.flush();
    }

    intmp.reset();
    
    infile.close();
    outfile.close();

    delete tree;
    return 0;
}
```
## Uncompress.cpp
```
#include <queue>
#include <vector>
#include <fstream>
#include "Helper.hpp"
#include "HCTree.hpp"
#include <algorithm>
using namespace std;

int main(int argc, char* argv[]){
    if(argc != 3){
        error("invalid number of arguments");
    }
    ifstream infile(argv[1]);
    ofstream outfile(argv[2]);
    //empty
    FancyInputStream intmp(argv[1]);
    if(intmp.filesize() == 0){
        return 0;
    }
    vector<int> freqs(256, 0);
    FancyInputStream fin(argv[1]);
    FancyOutputStream fout(argv[2]);
    //unoptimized header
    HCTree *tree = new HCTree();
    int sym = tree->helpdec(fin);
    int symtmp;
    int count = 0;
    unsigned char ch;

    while(1){
        symtmp = tree->helpdec(fin);
        count++;
        ch = (unsigned char)symtmp;
        outfile<<ch;
        //if all symbols decoded
        if(count == sym){
            break;
        }
    }
    //delete
    /*queue<HCNode*> nq;
    vector<HCNode*> nodes;
    HCNode *node;
    nq.push(root);
    while(1){
        node = nq.top();
        np.pop();
        nodes.push_back(node);
        if(node->c0 != nullptr){
            nq.push(node->c0);
        }
        if(node->c1 != nullptr){
            nq.push(node->c1);
        }
        if(nq.empty() == true){
            break;
        }
    }
    for(auto HCNode* n: nodes){
        delete n;
    } */

    infile.close();
    outfile.close();
    delete tree;
    return 0;
}
```
