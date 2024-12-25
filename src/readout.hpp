#ifndef READOUT_HPP
#define READOUT_HPP

#include <ostream>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <filesystem>
#include <fstream>
#include <algorithm> //for std::max
#include <unordered_map>
#include <mutex>

#include <fcntl.h>      // for open()
#include <unistd.h>     // for read() and close()
#include <sys/stat.h>   // for stat()
#include <cstring>      // for std::memcpy

#include "../../include/fast_float.h"
//#include "../../include/fast_convert.h"
//#include "../include/pdqsort.h"

//#include "../include/yas/mem_streams.hpp"
//#include "../include/yas/binary_iarchive.hpp"
//#include "../include/yas/binary_oarchive.hpp"
//#include "../include/yas/std_types.hpp"
//#include "../include/yas/file_streams.hpp"




struct star {
    std::string id;
    std::vector<double> x;
    std::vector<float> y;
    std::vector<float> dy;
    //float median;

    /*
    inline void estimate_median(){

    double y_sort[x.size()];
    for (int i=0; i<x.size(); ++i) {y_sort[i] = y[i];}
    pdqsort_branchless(&y_sort[0], &y_sort[x.size() - 1]);

    median = y_sort[int(x.size() * 0.5)];
    // std::cout << id << ".dat " << median << std::endl; //print approximate median value for file

    }
    */


inline void read(const std::string& in_file) {
    // Open file using POSIX open
    int fd = open(in_file.c_str(), O_RDONLY); // Open file in read-only mode
    if (fd == -1) { std::cerr << "File cannot be opened: " << strerror(errno) << std::endl; return;}

    // Get the size of the file
    struct stat file_stat;
    if (fstat(fd, &file_stat) == -1) {std::cerr << "Failed to get file size: " << strerror(errno) << std::endl; close(fd); return; }
    auto file_size = file_stat.st_size;

    // Allocate buffer of appropriate size
    char* dataBuffer = static_cast<char*>(malloc(file_size));
    if (dataBuffer == nullptr) { std::cerr << "Memory allocation failed" << std::endl; close(fd); return; }

    // Read the entire file into the buffer
    ssize_t bytes_read = ::read(fd, dataBuffer, file_size);
    if (bytes_read == -1) { std::cerr << "Failed to read file: " << strerror(errno) << std::endl; free(dataBuffer); close(fd); return;
    } else if (bytes_read == 0) { std::cout << "No data read from file." << std::endl; free(dataBuffer); close(fd); return;}

    size_t lines = 0;

    // Count lines in the loaded data
    for (size_t i = 0; i < bytes_read; ++i) {if (dataBuffer[i] == '\n') {++lines;}}

    // Set up iterators for parsing from dataBuffer
    const char* it = dataBuffer;
    const char* end = it + bytes_read;

    double tempX; float tempY; float tempDY;

    while (it < end) {
        auto temp = fast_float::from_chars(it, end, tempX);
        if (temp.ec != std::errc()) { std::cout << "Parsing failed" << std::endl; break; }
        if (temp.ptr[0] != ' ' && temp.ptr[0] != '\n') { std::cout << "Unexpected delimiter" << std::endl; }

        temp = fast_float::from_chars(temp.ptr + 1, end, tempY);
        if (temp.ec != std::errc()) { std::cout << "Parsing failed" << std::endl; break; }
        if (temp.ptr[0] != ' ' && temp.ptr[0] != '\n') { std::cout << "Unexpected delimiter" << std::endl; }

        temp = fast_float::from_chars(temp.ptr + 1, end, tempDY);
        if (temp.ec != std::errc()) { std::cout << "Parsing failed" << std::endl; break; }
        if (temp.ptr[0] != ' ' && temp.ptr[0] != '\n') { std::cout << "Unexpected delimiter" << std::endl; }

        it = temp.ptr + 1;

        x.push_back(tempX); y.push_back(tempY); dy.push_back(tempDY);
    }

    // Extract star id from filename
    id = in_file.substr(in_file.find_last_of('/') + 1, in_file.find_last_of('.') - in_file.find_last_of('/') - 1);

    // Free allocated memory
    free(dataBuffer);

    // Close the file descriptor
    close(fd);
}

/*

    template <typename Archive>
    inline void serialize(Archive& ar) {
    ar & id & x & y & dy;
    }

*/

/*
inline void read(const std::string& in_file) { //fast_convert method
    // Open file using POSIX open
    int fd = open(in_file.c_str(), O_RDONLY); // Open file in read-only mode
    if (fd == -1) {std::cerr << "File cannot be opened: " << strerror(errno) << std::endl; return;}

    // Get the size of the file
    struct stat file_stat;
    if (fstat(fd, &file_stat) == -1) {std::cerr << "Failed to get file size: " << strerror(errno) << std::endl; close(fd); return;}
    auto file_size = file_stat.st_size;

    // Allocate buffer of appropriate size using malloc
    char* dataBuffer = static_cast<char*>(malloc(file_size + 1)); // +1 for null-termination
    if (!dataBuffer) { std::cerr << "Memory allocation failed" << std::endl; close(fd); return;}

    // Read the entire file into the buffer
    ssize_t bytes_read = ::read(fd, dataBuffer, file_size);
    if (bytes_read == -1) {std::cerr << "Failed to read file: " << strerror(errno) << std::endl; free(dataBuffer); close(fd); return;
    } else if (bytes_read == 0) {std::cout << "No data read from file." << std::endl; free(dataBuffer); close(fd); return;}

    dataBuffer[bytes_read] = '\0'; // Null-terminate the buffer for safe string handling

    char* it = dataBuffer;
    char* end = dataBuffer + bytes_read;

    double tempX; float tempY; float tempDY;

    while (it < end) {
        // Parse tempX using fast_convert's fast_strtod
        tempX = fast_strtod(it, &it);
        if (it == nullptr || it >= end) {std::cout << "Parsing failed for X" << std::endl; break;}
        if (*it != ' ' && *it != '\n') {std::cout << "Unexpected delimiter after X" << std::endl; break;}

        it++; // Move past the delimiter

        // Parse tempY using fast_convert's fast_strtof
        tempY = fast_strtof(it, &it);
        if (it == nullptr || it >= end) {std::cout << "Parsing failed for Y" << std::endl; break;}
        if (*it != ' ' && *it != '\n') {std::cout << "Unexpected delimiter after Y" << std::endl; break;}

        it++; // Move past the delimiter

        // Parse tempDY using fast_convert's fast_strtof
        tempDY = fast_strtof(it, &it);
        if (it == nullptr || it >= end) {std::cout << "Parsing failed for DY" << std::endl; break;}
        if (*it != ' ' && *it != '\n') {std::cout << "Unexpected delimiter after DY" << std::endl; break;}

        it++; // Move to the next token

        // Store the parsed values into vectors
        x.push_back(tempX); y.push_back(tempY); dy.push_back(tempDY);
    }

    // Free allocated memory
    free(dataBuffer);

    // Extract the star ID from the file name
    id = in_file.substr(in_file.find_last_of('/') + 1, in_file.find_last_of('.') - in_file.find_last_of('/') - 1);

    // Close the file descriptor
    close(fd);
}
*/



inline void linreg() {
    double tmp = 0;

    // Center the measurement times - increases precision of future computation
    for (unsigned int i = 0; i < x.size(); i++) {tmp += x[i];}  tmp /= double(x.size());
    for (unsigned int i = 0; i < x.size(); i++) {x[i] -= tmp;}

    // Initialize sums for weighted regression
    double sumx = 0, sumxsq = 0, sumy = 0, sumxy = 0, sumw = 0; // Add sum of weights
    double lin = 0, c = 0, w;

    for (unsigned int i = 0; i < x.size(); i++) {
        w = 1 / (dy[i] * dy[i]); // weight as the inverse of variance
        if (std::isnormal(w) && w > 0){
            sumw += w;               // accumulate total weights
            sumx += x[i] * w;       // weighted sum of x
            sumxsq += x[i] * x[i] * w; // weighted sum of x^2
            sumy += y[i] * w;       // weighted sum of y
            sumxy += (x[i] * y[i]) * w; // weighted sum of x*y
        }
        else {dy[i] = 999.9;}
    }

    // Calculate the denominator for the slope and intercept
    double denum = (sumw * sumxsq) - (sumx * sumx);

    // Compute the slope (lin) and intercept (c)
    lin = ((sumw * sumxy) - (sumx * sumy)) / denum;
    c = ((sumy * sumxsq) - (sumx * sumxy)) / denum;

    // Adjust the y values based on the regression line
    for (unsigned int i = 0; i < y.size(); i++) {y[i] -= (lin * x[i]) + c;}
    }

};


/*
struct photometry {
    std::vector<star> stars;
    std::unordered_map<std::string, int> id;
    std::mutex mutex;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & stars & id;
    }


    inline void add(const star &star) {
    std::lock_guard<std::mutex> lock(mutex);

    stars.push_back(star);
    id[star.id] = stars.size() - 1;
    }

    // Deserialize and load the photometry struct from file
    void load(const std::string& file) {
        const char* file_cstr = file.c_str();

        yas::file_istream fis(file_cstr);
        yas::binary_iarchive<yas::file_istream> ia(fis);

        ia & *this;
    }

    void load_dat(std::string in_dir){
        std::vector<std::filesystem::path> files;
        for (auto& entry : std::filesystem::directory_iterator(in_dir)) {
            if (entry.is_regular_file() && entry.path().extension() == ".dat") {
                files.push_back(entry.path());
            }
        }

        #pragma omp parallel for
        for (unsigned int i = 0; i < files.size(); i++) {
            star data;
            data.read(files[i]);
            add(data);}
    }
};
*/


#endif
