/* Created by the Pantazes Lab at Auburn University.
 *
 * This file implements the general list directory method */

// This file is supposed to be included by Methods.h
#ifndef Methods_Loading_Status
#error Methods::listdir.h must be included by Methods.h
#endif

// Implement the method
// this is the function to list the files in a directory
vector<string> METHODS::listdir(string dir) {
    vector<string> files;
    string cmd = "ls " + dir;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        cout<<"Error: could not open pipe"<<endl;
        // exit the program
        exit(1);
    }
    char buffer[128];
    while (!feof(pipe)) {
        if (fgets(buffer, 128, pipe) != NULL) {
            files.push_back(string(buffer));
        }
    }
    pclose(pipe);
    // remove the newline characters
    for (size_t i=0; i<files.size(); ++i) {
        files[i].erase(files[i].size() - 1);
    }
    return files;
}              
