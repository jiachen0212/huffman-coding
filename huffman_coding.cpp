/**
 * This is an implementation of Huffman coding.
 *
 * The core algorithm is taken from the CLR book (Introduction of Algorithms),
 * Chapter 16.3, and directly used to implement the 'build_tree()' routine.
 *
 * After the tree is built, a code table that maps a character to a binary
 * code is built from the tree, and used for encoding text. Decoding is done
 * by traversing the Huffman tree, as prescribed by the algorithm.
 *
 * Binary codes are represented by std::vector<bool>, which is a specialized
 * vector that optimizes space.
 */

#include <vector>
#include <sstream>
#include <queue>
#include <map>
#include <string.h>
#include <algorithm>
#include <string>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <time.h>
#include <fstream>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>
#include <typeinfo>
#include <thread>
#include <future>      //std::future std::promise
#include <utility>     //std::ref
#include <chrono>
#include <unistd.h>
using namespace std;


// A Huffman Tree Node
struct HuffmanTree {
    char c; // character in an alphabet
    int cfreq; // frequency of c.
    struct HuffmanTree *left;
    struct HuffmanTree *right;
    HuffmanTree(char c, int cfreq, struct HuffmanTree *left=NULL,
                struct HuffmanTree *right=NULL) :
        c(c), cfreq(cfreq), left(left), right(right) {
    }
    ~HuffmanTree() {
        delete left, delete right;
    }
    // Compare two tree nodes
    class Compare {
    public:
        bool operator()(HuffmanTree *a, HuffmanTree *b) {
            return a->cfreq > b->cfreq;
        }
    };
};

/**
 * Builds a Huffman Tree from an input of alphabet C, where C is a vector
 * of (character, frequency) pairs.
 */
HuffmanTree *build_tree(vector< pair<char, unsigned> > &alph) {
    // First build a min-heap
    // Build leaf nodes first
    priority_queue<HuffmanTree *, vector<HuffmanTree *>, HuffmanTree::Compare > alph_heap;
    for (vector< pair<char, unsigned> >::iterator it = alph.begin();
         it != alph.end(); ++it) {
        HuffmanTree *leaf = new HuffmanTree(it->first, it->second);
        alph_heap.push(leaf);
    }

    // HuffmanTree algorithm: Merge two lowest weight leaf nodes until
    // only one node is left (root).
    HuffmanTree *root = NULL;
    while (alph_heap.size() > 1) {
        HuffmanTree *l, *r;
        l = alph_heap.top();
        alph_heap.pop();
        r = alph_heap.top();
        alph_heap.pop();
        root = new HuffmanTree(0, l->cfreq + r->cfreq, l, r);
        alph_heap.push(root);
    }

    return root;
}

/**
 * Prints the tree nodes in breadth-first order
 */
void print_tree(HuffmanTree *t) {
    deque< pair<HuffmanTree *, int> > q;

    q.push_back(make_pair(t, 0));
    int curlevel = -1;
    while (!q.empty()) {
        HuffmanTree *parent = q.front().first;
        int level = q.front().second;
        q.pop_front();
        if (curlevel != level) {
            curlevel = level;
            cout << "Level " << curlevel << endl;
        }
        cout << parent->cfreq << " " << parent->c << endl;
        if (parent->left)
            q.push_back(make_pair(parent->left, level + 1));
        if (parent->right)
            q.push_back(make_pair(parent->right, level + 1));
    }
}

typedef vector<bool> code_t;
typedef map<char, code_t> codetable;

/**
 * Makes a lookup table (std::map) of (c -> code) from a HuffmanTree, where
 * code is an unsigned long representing the binary code.
 */
map<char, code_t> build_lookup_table(HuffmanTree *htree) {
    codetable lookup;  // 二进制表  char -> code
    deque< pair<HuffmanTree *, code_t> > q;

    q.push_back(make_pair(htree, code_t()));
    while (!q.empty()) {
        HuffmanTree *node, *lc, *rc;
        code_t code;
        node = q.front().first;
        code = q.front().second;
        q.pop_front();
        lc = node->left;
        rc = node->right;
        if (lc) {
            // HuffmanTree is always full (either no children or two children)
            // Left child is appended a 0 and right child a 1.
            code_t code_cp(code);
            q.push_back(make_pair(lc, (code.push_back(0), code)));
            q.push_back(make_pair(rc, (code_cp.push_back(1), code_cp)));
        } else {
            // Leaf node: contains the character
            lookup.insert(make_pair(node->c, code));
            cout << "(" << node->c << ", ";
            for (unsigned i = 0; i < code.size(); i++) {
                cout << code[i];
            }
            cout << ")" << endl;
        }
    }

    return lookup;
}

/**
 * Encodes an input string. returns a byte vector.
 */
code_t encode(string input, codetable &lookup, code_t &result) {

    // code_t reslut
    string::iterator it = input.begin();
    string::iterator it_end = input.end();
    for (it; it != it_end; ++it) {
            code_t b = lookup[*it];
            result.insert(result.end(), b.begin(), b.end());
    }
    // return result;   // 编码得到的01 byte 序列 vector存储...
}

/**
 * Look up the next valid code in @biter using @htree and returns the
 * resulting string. Note the iterator @biter is advanced by the actual
 * length of the next valid code, which varies.
 */

char code_lookup(code_t::iterator &biter, const code_t::iterator &biter_end,
                 const HuffmanTree *htree) {
    const HuffmanTree *node = htree;

    while (true) {
        if (!node->left)
        {   // Huffman tree is full: always contains both children or none.
            break;
        }

        if (biter == biter_end) {
            throw std::out_of_range("No more bits");
        }

        if (*biter) {
            node = node->right;
        } else {
            node =node->left;
        }
        ++biter;
    }

    return node->c;
}

/**
 * Decodes a compressed string represented by a bit vector (vector<char>)
 * @compressed, using a HuffmanTree @htree.
 * Returns the original string.
 */
string decode(code_t &compressed, const HuffmanTree *htree)
{
    string result;
    code_t::iterator biter = compressed.begin();
    code_t::iterator biter_end = compressed.end();
    while (true) {
        try {
            result += code_lookup(biter, biter_end, htree);
        } catch (const std::out_of_range &oor) {
            // Iterator exhausted.
            break;
        }
    }
    return result;
}

/**
 * Tests
 */
// Make frequency table from a string.
vector< pair<char, unsigned> > make_freq_table(string inp) {
    map<char, unsigned> cfmap;
    vector< pair<char, unsigned> >cfvec;

    for (unsigned i = 0; i < inp.size(); i++) {
        if (cfmap.find(inp[i]) == cfmap.end()) {
            cfmap.insert(make_pair(inp[i], 1));
        }
        cfmap[inp[i]] += 1;
    }

    for (map<char, unsigned>::iterator it = cfmap.begin();
         it != cfmap.end(); ++it) {
        cfvec.push_back(make_pair(it->first, it->second));
    }

    return cfvec;
}

string bitvec_to_string(code_t &bitvec) {
    string result;
    size_t nbits;

    nbits = bitvec.size() & 7;  //4

    // Write the number of "hanging bits" at the first byte
    result += static_cast<char>(nbits); // at most 7

    char byte = 0;
    int bitvec_size = bitvec.size();
    for (unsigned i = 0; i < bitvec_size; i++) {
        unsigned boff = i & 7;
        byte |= bitvec[i] << boff;
        if (boff == 7) {
            // Write a byte
            result += byte;
            byte = 0;
        }
    }

    if (nbits) {
        result += byte;
    }
    // cout << bitvec_size << "----" << result.size() << endl;
    return result;
}

code_t string_to_bitvec(string packed) {
    code_t result;
    int packed_size = packed.size();
    // assert(packed.size());
    if (packed_size == 1) {
        assert(packed[0] == 0);
        return result;
    }
    unsigned nbits = packed[0];
    string::iterator it = packed.begin() + 1;
    string::iterator packed_end = packed.end();
    for (it; it != packed_end; ++it) {
        for (unsigned i = 0; i < 8; i++) {
            result.push_back((*it >> i) & 1);
        }
    }
    // fix the last byte
    if (nbits) {
        for (unsigned i = 0; i < (8 - nbits); i++) {
            result.pop_back();
        }
    }

    return result;
}

void encoder(string &s, codetable &ctbl, code_t &t)
{
    // Encode
    // timeval start, end;
    // float time_use=0;
    // gettimeofday(&start, NULL);
    encode(s, ctbl, t);
    // gettimeofday(&end, NULL);
    // time_use=(end.tv_sec-start.tv_sec)*1000000 +(end.tv_usec-start.tv_usec);
    // printf("Single thread encode cost: %f秒\n",time_use/1000000);
    // cout << "encoded (compression ratio: " << 1 - ((float)(t.size() / 8) / s.size()) << ")" << endl;
}

void decoder(string &ans, string &c, const HuffmanTree *htree)
{
    // Decode
    // timeval start1, end1;
    // float time_use1=0;
    // gettimeofday(&start1, NULL);
    // string packed = bitvec_to_string(t);
    // cout << t.size() << " " << packed.size() << endl;
    code_t t1 = string_to_bitvec(c);
    // cout << t.size() << " " << t1.size() << endl;
    // assert(std::equal(t.begin(), t.end(), t1.begin()));
    ans = decode(t1, htree);  // 这句有点坑  string ans = decode(t1, htree) 报奇怪内存错误...
    // gettimeofday(&end1, NULL);
    // time_use1=(end1.tv_sec-start1.tv_sec)*1000000 +(end1.tv_usec-start1.tv_usec);
    // printf("Single thread, decode cost: %f秒\n",time_use1/1000000);
    //// assert(s1 == s);
    //// delete htree;
    // cout << "string size:" << c.size() << "\tdecompressed size: " << ans.size() << endl;
}

string readFileIntoString(){
    ifstream ifile("compressed.txt");
    ostringstream buf;
    char ch;
    while(buf&&ifile.get(ch))
    buf.put(ch);
    return buf.str();
}

int main() {
    string filename("y.txt");
    string text;
    fstream in(filename.c_str());
    string s;
    while(in>>text)
    {
        s += text;
    }
    int n_thread = 12;
    int sub_s = s.size() / n_thread;
    vector<string> sub_strings;
    for (unsigned i=0; i<n_thread; i++ ){
        sub_strings.push_back(s.substr(i*sub_s, sub_s));
    }
    // code_t t_1, t_2, t_3, t_4, t_5, t_6, t_7, t_8, t_9, t_10;   // 存储编码各子线程结果
    vector<code_t> t_s(n_thread);

    // 编解码总计时:
    timeval start, end;
    float time_use=0;
    gettimeofday(&start, NULL);

    // 构建 dict huffman tree
    vector< pair<char, unsigned> > cfvec = make_freq_table(s);
    HuffmanTree *htree = build_tree(cfvec);    //print_tree(htree);
    codetable ctbl = build_lookup_table(htree);

    // encoder 并发
    thread threads[n_thread];
    for (unsigned i=0; i<n_thread; i++){
        threads[i] = thread(encoder, ref(sub_strings[i]), ref(ctbl), ref(t_s[i]));
    }
    for (auto &thread : threads)
        thread.join();

    // save compressed file
    ofstream file("compressed.txt");
    vector<int> ss_start_index;  // 从compressed.txt中截取各个子编码string,每个子模块的起点index存起来...
    vector<int> ss_lens;
    string tmp;
    int tmp_size;
    int start_index = 0;
    // for (unsigned i=0; i<n_thread; i++){
    for (auto &ti : t_s){
        // tmp = bitvec_to_string(t_s[i]);
        tmp = bitvec_to_string(ti);
        tmp_size = tmp.size();
        ss_lens.push_back(tmp_size);
        ss_start_index.push_back(start_index);
        start_index += (tmp_size+1);  // 加1是因为每个string后面默认的'\0'
        file << tmp << endl;  // 只能依次写入.. 因为不同sub_thread除以8的余数与总code_t除以8的余数，合起来是有差别的...
    }

    /************   Decode   ************/
    
    // to get sub compressed size
    int all = ss_start_index[n_thread-1] + ss_lens[n_thread-1] + 1;
    // load compressed txt file
    string c; // string 后面的那个'\0'很fan ren
    c=readFileIntoString();

    vector<string> c_s;
    vector<string> ress(n_thread);

    // 解码计时:
    timeval start1, end1;
    float time_use1=0;
    gettimeofday(&start1, NULL);

    for (unsigned i=0; i<n_thread; i++){
        c_s.push_back(c.substr(ss_start_index[i], ss_lens[i]));
    }
    // Decode 并发:
    thread de_threads[n_thread];
    for (unsigned i=0; i<n_thread; i++){
        de_threads[i] = thread(decoder, ref(ress[i]), ref(c_s[i]), htree);
    }
    for (auto &thread : de_threads)
        thread.join();

    // concate results
    // string s_res = res1 + res2 + res3 + res4 + res5 + res6 + res7 + res8 + res9 + res10;  // 解码结果
    string s_res;
    for (auto &res_ : ress)
        s_res += res_;

    // only decompress time
    gettimeofday(&end1, NULL);
    time_use1=(end1.tv_sec-start1.tv_sec)*1000000 +(end1.tv_usec-start1.tv_usec);
    printf("解码耗时: %f秒\n",time_use1/1000000);

    // en+decode time
    gettimeofday(&end, NULL);
    time_use=(end.tv_sec-start.tv_sec)*1000000 +(end.tv_usec-start.tv_usec);
    printf("编解码总耗时: %f秒\n",time_use/1000000);
    cout << "compression ratio: " << 1 - ((float) all / s.size()) << endl;
    printf("compreassed file %f kb\n",(float) all / 1024);

    // 解码恢复数据 check
    assert(s_res == s);
    delete htree;
}
