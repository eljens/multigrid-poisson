#include "../include/parser.h"

Settings parser(int argc, char * argv[]){
    string str;
    Settings set;
    for (int i=1;i<argc;i++){
        str = argv[i];
        if (str.compare("-x") == 0){
            set.dims[0] = stoi(argv[i+1]);
            i++;
        }
        else if (str.compare("-y") == 0){
            set.dims[1] = stoi(argv[i+1]);
            i++;
        }
        else if (str.compare("-z") == 0){
            set.dims[2] = stoi(argv[i+1]);
            i++;
        }
        else if ((str.compare("-tolerance") == 0) || (str.compare("-tol") == 0)){
            set.tolerance = stod(argv[i+1]);
            i++;
        }
        else if (str.compare("-maxiter") == 0){
            set.maxiter = stoi(argv[i+1]);
            i++;
        }
        else if ((str.compare("-levels") == 0) || (str.compare("-l") == 0)){
            set.levels = stoi(argv[i+1]);
            i++;
        }
        else {
            cerr << "Unknown input " << str << endl;
        }
    }
    return set;
}
