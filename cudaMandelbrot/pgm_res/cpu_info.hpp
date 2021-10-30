#include<iostream>
#include<string>
#include<fstream>

std::string get_cpu_info( std::string pattern = "model name"){
    std::ifstream cpu_info;
    std::string line;
    cpu_info.open("/proc/cpuinfo");
    while(getline(cpu_info,line)){
        if(line.substr(0,pattern.size()) == pattern){
            //std::cout << line << std::endl;
            break;
        }
    }

    size_t pos{0};
    for(int i = 0; i<line.size(); ++i){
        if(line[i] == ':'){
            pos = i;
            break;
        }
    }
    cpu_info.close();
    //std::cout << line.substr(pos + 2,line.size()) << std::endl;
    return line.substr(pos + 2,line.size());
    
}