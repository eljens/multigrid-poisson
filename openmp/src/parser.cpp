#include "parser.h"

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
        else if ((str.compare("-save") == 0) || (str.compare("-s") == 0)){
            set.print_result = stoi(argv[i+1]);
            i++;
        }
        else if ((str.compare("-length") == 0)){
            set.lengthx = stod(argv[i+1]);
            i++;
        }
        else if ((str.compare("-stats") == 0)){
            set.stats_file = argv[i+1];
            set.write_final_stats = true;
            i++;
        }
        else {
            cerr << "Unknown input " << str << endl;
        }
    }
    set.h = set.lengthx/((double_t) set.dims[0] - 1.0);
    return set;
}
