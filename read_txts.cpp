#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

/*
读入需编码的test.txt数据.
test.txt: -10,100, 255, ...-266
*/

int main()
{
    ifstream inf;
    inf.open("test.txt", ifstream::in);
    string line;
    size_t first_index = 0;
    size_t second_index = 0;
    vector<int> cc;  // 接收 .txt 中的数据流..

    while (!inf.eof())
    {
        getline(inf,line);
        int len_ = line.size();
        // cout << line << endl;
        first_index = line.find(',',0);  // 从 index 0 开始查找, 且找到的第一个','的index是first_index..
        // cout << "first_index " << first_index << endl;
        int num = atoi(line.substr(0,first_index).c_str());  // atoi() string转int
        cc.push_back(num);
        while (first_index < len_)
        {
            second_index = line.find(',', first_index + 1);
            // cout << "second_index " << second_index << endl;
            int num_ = atoi(line.substr(first_index + 1,second_index-first_index-1).c_str());
            cc.push_back(num_);
            first_index = second_index;
        }
    }
    cc.pop_back();  // 删除最末字符,它不是txt中的元素

    inf.close();
    for (auto &c : cc)
        cout << c << endl;
    return 0;
}