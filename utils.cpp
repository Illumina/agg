#include "utils.h"
using namespace std;

int *zeros(int n) {
    int * ret = new int[n];
    for(int i = 0; i < n; i++) ret[i] = 0;
    return(ret);
}

bool fileexists(const string& fname){
    ifstream ifile(fname.c_str());
    return ifile;
}

string percent(int num, int den) {
    char buffer[100];
    sprintf(buffer, "%d / %d\t( %.4f%% )\t", num, den, 100. * float(num) / float(den));
    return(buffer);
}

/*void die(const string& s) {
    cerr << "ERROR: " << s << "\nExiting..." << endl;
    exit(1);
}*/

int strsplit(const string& input, const char split, vector<string> & out) {
    istringstream ss(input);
    out.clear();
    string tmp;
    int count=0;
    while(std::getline(ss, tmp, split))  {
        out.push_back(tmp);    
        count++;
    }
    return(count);
}

string join(const vector<string>& input, const string& delim) {
    string ret = input[0];
    for(int i=1;i<input.size();i++) {
        ret += ",";
        ret += input[i];
    }
    return(ret);
}

vector<int> match(const vector<string>& x, const vector<string>& y) {
    vector<int> ret;
    for(int i = 0; i < y.size(); i++) {
        int j = 0;
        while(x[j] != y[i]) {
            j++;
            if(j >= x.size()) {
                die("sample in y is not in x");
            }
        }
        ret.push_back(j);
    }
    return(ret);
}

/*void warn(const string& s)  {
    cerr << "WARNING: "<<s<<endl;
}*/
